################################################################################
## Script to run:
##  - Lifestyle volcano (Cohen's d, Bonferroni)
##  - Lifestyle × Sex interaction scan + plots
## using Data_for_NMR_ANALYSIS + some other metabolite_cols
################################################################################

## ─────────────────────────── 0. Setup ─────────────────────────────────────────
setwd("/Users/bm0211/RegEx/Water_Energy/2025_Metabolomic_Analysis/NMR_metabolomics/")
## Packages needed
library(rio)        # import()
library(dplyr)      # select, mutate, filter, %>%
library(tidyr)      # drop_na
library(purrr)      # map_dfr
library(tibble)     # tibble()
library(ggplot2)    # plotting
library(ggrepel)    # geom_text_repel
library(emmeans)    # emmeans()
library(ggpubr)     # ggarrange()
library(glue)       # glue() in error messages

################################################################################
## 1. Load data containing the Framingham score and NMR data and build Data_for_NMR_ANALYSIS
################################################################################

# This is the file produced after running your Python Framingham script
RAW_Data_for_NMR_ANALYSIS <- import(
  "/Users/bm0211/RegEx/Water_Energy/2025_Metabolomic_Analysis/RAW_Data_for_NMR_ANALYSIS.csv")

# Keep only rows with valid lifestyle info
RAW_Data_for_NMR_ANALYSIS <- RAW_Data_for_NMR_ANALYSIS %>%
  filter(!is.na(MW_scores_lifestyle), MW_scores_lifestyle != "")

# Collapse lifestyle: Pastoralist + PeriUrban = Non_Urban, Urban = Urban and LEVEL
RAW_Data_for_NMR_ANALYSIS <- RAW_Data_for_NMR_ANALYSIS %>%
  mutate(lifestyle_group_MW = ifelse(
      MW_scores_lifestyle %in% c("Pastoralist", "PeriUrban"),
      "Non_Urban", "Urban"))

RAW_Data_for_NMR_ANALYSIS$lifestyle_group_MW <- factor(
  RAW_Data_for_NMR_ANALYSIS$lifestyle_group_MW,
  levels = c("Non_Urban", "Urban"))

################################################################################
## 2. Subset to NMR + key covariates and preprocess
################################################################################
colnames(RAW_Data_for_NMR_ANALYSIS)
# Same column selection you used in the main script
columns_to_keep <- c(1,25,31,590, 49,600,602,875, 888,890,891, 892, 903, 905, 625:874, 589,489:509)

RAW_Data_for_NMR_ANALYSIS <- RAW_Data_for_NMR_ANALYSIS %>% dplyr::select(all_of(columns_to_keep))

# Remove Nightingale prefix from NMR variables
colnames(RAW_Data_for_NMR_ANALYSIS) <- gsub("^Nightingale_NMR_", "", colnames(RAW_Data_for_NMR_ANALYSIS))

# Drop percentage columns
RAW_Data_for_NMR_ANALYSIS <- RAW_Data_for_NMR_ANALYSIS %>% dplyr::select(-matches("_pct$"), -matches("(?i)\\.PCT$"))

# Numeric columns to log+scale (MAKE SURE same indices as in your DATA)
columns_to_process <- c(2, 8:11, 14:193)
# Ensure covariates are in the right format
RAW_Data_for_NMR_ANALYSIS <- RAW_Data_for_NMR_ANALYSIS %>%
  mutate(lifestyle_group_MW = recode(MW_scores_lifestyle,
                                     "Pastoralist" = "Non_Urban",
                                     "PeriUrban"   = "Non_Urban",
                                     "Urban"       = "Urban") |> 
           factor(levels = c("Non_Urban", "Urban")))

RAW_Data_for_NMR_ANALYSIS$Sex <- as.factor(RAW_Data_for_NMR_ANALYSIS$Sex)
RAW_Data_for_NMR_ANALYSIS$lifestyle_group_MW <- as.factor(
  RAW_Data_for_NMR_ANALYSIS$lifestyle_group_MW)

## Helper: convert selected columns to numeric
convert_to_numeric <- function(df, columns) {
  for (col in columns) {
    nm <- colnames(df)[col]
    df[[nm]] <- as.numeric(as.character(df[[nm]]))
  }
  df}

RAW_Data_for_NMR_ANALYSIS <- convert_to_numeric(RAW_Data_for_NMR_ANALYSIS, columns_to_process)

## Helper: remove outliers via IQR rule
remove_outliers_iqr <- function(df, columns) {
  for (col in columns) {
    nm <- colnames(df)[col]
    if (is.numeric(df[[nm]])) {
      Q1 <- quantile(df[[nm]], 0.25, na.rm = TRUE)
      Q3 <- quantile(df[[nm]], 0.75, na.rm = TRUE)
      IQR <- Q3 - Q1
      lower <- Q1 - 1.5 * IQR
      upper <- Q3 + 1.5 * IQR
      df[[nm]][df[[nm]] < lower | df[[nm]] > upper] <- NA
    }
  }
  df
}

# Outlier removal in raw space
RAW_Data_for_NMR_ANALYSIS <- remove_outliers_iqr(
  RAW_Data_for_NMR_ANALYSIS, columns_to_process)

## Log10-transform and z-score
standardize_columns <- function(df, columns) {
  for (col in columns) {
    nm <- colnames(df)[col]
    log_values <- log10(df[[nm]])
    df[[nm]] <- scale(log_values, center = TRUE, scale = TRUE)
  }
  df
}

Data_for_NMR_ANALYSIS <- standardize_columns(RAW_Data_for_NMR_ANALYSIS, columns_to_process)

# Outlier removal on transformed data
Data_for_NMR_ANALYSIS <- remove_outliers_iqr(Data_for_NMR_ANALYSIS, columns_to_process)

# Make sure Age is numeric and Sex is factor in the final analysis object
Data_for_NMR_ANALYSIS$Age <- as.numeric(Data_for_NMR_ANALYSIS$Age)
Data_for_NMR_ANALYSIS$Sex <- factor(Data_for_NMR_ANALYSIS$Sex)

################################################################################
## Define metabolite columns
################################################################################

# These are the metabolites you scan for lifestyle/sex/age effects
metabolite_cols <- names(Data_for_NMR_ANALYSIS[, c(9:13, 15:56, 58:64, 66:68, 70:82, 84:86, 89:189)])

################################################################################
## analysis block
################################################################################

## lock lifestyle order again just to make sure (Level1 = Non_Urban, Level2 = Urban)
Data_for_NMR_ANALYSIS$lifestyle_group_MW <- factor(
  Data_for_NMR_ANALYSIS$lifestyle_group_MW, levels = c("Non_Urban", "Urban"))

## ------------------------ LOOP: fit + Cohen's d ------------------------
res_list <- vector("list", length(metabolite_cols))
names(res_list) <- metabolite_cols
for (m in metabolite_cols) {
  y    <- Data_for_NMR_ANALYSIS[[m]]
  Age  <- Data_for_NMR_ANALYSIS$Age
  Sex  <- Data_for_NMR_ANALYSIS$Sex
  LSG  <- Data_for_NMR_ANALYSIS$lifestyle_group_MW
  
  ok <- complete.cases(y, Age, Sex, LSG)
  if (!any(ok)) next
  
  y   <- y[ok]; Age <- Age[ok]; Sex <- Sex[ok]; LSG <- LSG[ok, drop=TRUE]
  
  ## need both groups and enough data
  if (length(unique(LSG)) < 2 || length(y) < 200) next
  if (any(table(LSG) < 2))  next  # need variance in each group
  
  ## fit model: y ~ lifestyle + Sex + Age
  fit <- lm(y ~ LSG + Sex + Age)
  co  <- summary(fit)$coefficients
  ## grab the Urban contrast row
  beta <- se <- pval <- NA_real_
  if ("LSGUrban" %in% rownames(co)) {
    beta <- co["LSGUrban", "Estimate"]
    se   <- co["LSGUrban", "Std. Error"]
    pval <- co["LSGUrban", "Pr(>|t|)"]}
  
  ## Cohen's d = (mean(Urban) - mean(Non_Urban)) / pooled SD
  g1 <- y[LSG == "Non_Urban"]; g2 <- y[LSG == "Urban"]
  n1 <- length(g1); n2 <- length(g2)
  d  <- NA_real_
  if (n1 > 1 && n2 > 1 && sd(g1) > 0 && sd(g2) > 0) {
    sp <- sqrt(((n1 - 1)*var(g1) + (n2 - 1)*var(g2)) / (n1 + n2 - 2))
    d  <- (mean(g2) - mean(g1)) / sp  # Urban − Non_Urban
  }
  
  res_list[[m]] <- data.frame(
    Metabolite        = m,
    term              = "lifestyle_group_MWUrban",
    estimate          = beta,
    std.error         = se,
    p.value           = pval,
    N                 = length(y),
    n_Non_Urban       = n1,
    n_Urban           = n2,
    mean_Non_Urban    = if (n1) mean(g1) else NA_real_,
    mean_Urban        = if (n2) mean(g2) else NA_real_,
    sd_Non_Urban      = if (n1 > 1) sd(g1) else NA_real_,
    sd_Urban          = if (n2 > 1) sd(g2) else NA_real_,
    EffectSize_d      = d,
    Level1            = "Non_Urban",
    Level2            = "Urban",
    stringsAsFactors  = FALSE
  )
}

## bind
all_coefficients <- do.call(rbind, res_list[!vapply(res_list, is.null, logical(1))])
if (is.null(all_coefficients) || nrow(all_coefficients) == 0L) {
  stop("No valid metabolites after filtering.")
}

## FDR + tidy volcano columns
all_coefficients$adjusted_p_values <- p.adjust(all_coefficients$p.value, method = "bonferroni")
all_coefficients$nlog10FDR <- -log10(all_coefficients$adjusted_p_values)
colnames(all_coefficients)
## keep only rows we plot (the Urban contrast + finite values)
volcano_df <- subset(
  all_coefficients,
  term == "lifestyle_group_MWUrban" & is.finite(EffectSize_d) & is.finite(nlog10FDR))

## --- build volcano classes + choose labels ---
sig <- volcano_df$adjusted_p_values < 0.05

volcano_df$volc_class <- ifelse(!sig, "NonSig",
                                ifelse(volcano_df$EffectSize_d > 0, "Urban_up", "NonUrban_up"))
volcano_df$volc_class <- factor(volcano_df$volc_class, levels = c("NonSig","Urban_up","NonUrban_up"))

# top each side among significant points
pos_idx <- which(sig & volcano_df$EffectSize_d > 0)
neg_idx <- which(sig & volcano_df$EffectSize_d < 0)
top_pos <- pos_idx[order(volcano_df$EffectSize_d[pos_idx], decreasing = TRUE)][seq_len(min(10, length(pos_idx)))]
top_neg <- neg_idx[order(volcano_df$EffectSize_d[neg_idx], decreasing = FALSE)][seq_len(min(10, length(neg_idx)))]
lab_idx <- c(top_pos, top_neg)

urban_col     <- "#FC8D62"
non_urban_col <- "green"

## --- ggplot volcano ---
ggplot(volcano_df, aes(x = EffectSize_d, y = nlog10FDR, color = volc_class)) +
  geom_point(alpha = 0.85, size = 2.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_text_repel(
    data = volcano_df[lab_idx, ],
    aes(label = Metabolite),
    size = 4,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2,
    min.segment.length = 0,
    segment.color = "grey55" ) +
  scale_color_manual(
    values = c(
      NonSig      = "grey80",
      Urban_up    = urban_col,      # significant & d > 0 (Urban higher)
      NonUrban_up = non_urban_col   # significant & d < 0 (Non_Urban higher)
    ), guide = "none") + 
  labs( title = "Volcano: Lifestyle contrast on Metabolites",
        x = "Cohen's d (Urban − Non_Urban)",
        y = "-log10(Bonferroni p)") + theme_minimal(base_size = 22)

## ------------------------ EXPORT ------------------------
write.csv(all_coefficients,
          "NMR_all_coefficients_results_SIMPLE.csv",
          row.names = FALSE)


#####################################################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################
# Checking for interaction effects
table(Data_for_NMR_ANALYSIS$lifestyle_group_MW)
# --- Evaluate models and store results ---
# Evaluate models and store results, skipping problematic ones
lrt_summary <- map_dfr(metabolite_cols, function(metab) {
  tryCatch({
    df <- Data_for_NMR_ANALYSIS %>%
      dplyr::select(Metabolite = !!sym(metab), lifestyle_group_MW, Sex, Age) %>% drop_na()
    
    # Check factor levels
    if (n_distinct(df$lifestyle_group_MW) < 2 | n_distinct(df$Sex) < 2) {
      stop("Not enough levels in factors for contrasts")
    }
    full_model <- lm(Metabolite ~ lifestyle_group_MW * Sex + Age, data = df)
    reduced_model <- lm(Metabolite ~ lifestyle_group_MW + Sex + Age, data = df)
    lrt <- anova(reduced_model, full_model)
    p_val <- lrt$`Pr(>F)`[2]
    tibble(Metabolite = metab, LRT_p = round(p_val, 6),
           Better_Model = ifelse(p_val < 0.05, "Full (with interaction)", "Reduced (no interaction)")
    )
  }, error = function(e) {
    message(glue::glue("Skipping {metab} due to error: {e$message}"))
    NULL
  })
})

# --- Count how many times each model was better ---
full_better <- sum(lrt_summary$Better_Model == "Full (with interaction)")
reduced_better <- sum(lrt_summary$Better_Model == "Reduced (no interaction)")
total <- full_better + reduced_better
cat("\nSummary:\n")
cat(paste("🧪 Full model better in", full_better, "metabolites (",
          round(100 * full_better / total, 1), "% )\n"))
cat(paste("⚪ Reduced model better in", reduced_better, "metabolites (",
          round(100 * reduced_better / total, 1), "% )\n"))
# --- Plot interaction effects for those where full model was better ---
int_metabs <- lrt_summary %>%
  filter(Better_Model == "Full (with interaction)") %>% pull(Metabolite)
em_all <- map_dfr(int_metabs, function(metab) {
  df <- Data_for_NMR_ANALYSIS %>%
    dplyr::select(Metabolite = !!sym(metab), lifestyle_group_MW, Sex, Age) %>% drop_na()
  fit <- lm(Metabolite ~ lifestyle_group_MW * Sex + Age, data = df)
  em <- emmeans(fit, ~ lifestyle_group_MW * Sex) %>% as.data.frame()
  em %>% mutate(Metabolite = metab)})

# Generate list of ggplots for significant interactions
interaction_plots <- map(int_metabs, function(metab) {
  df <- Data_for_NMR_ANALYSIS %>%
    dplyr::select(Metabolite = !!sym(metab), lifestyle_group_MW, Sex, Age) %>% drop_na()
  fit <- lm(Metabolite ~ lifestyle_group_MW * Sex + Age, data = df)
  em <- emmeans(fit, ~ lifestyle_group_MW * Sex) %>% as.data.frame()
  ggplot(em, aes(x = lifestyle_group_MW, y = emmean, colour = Sex, group = Sex)) +
    geom_point(position = position_dodge(0.2), size = 2) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                  width = .1, position = position_dodge(0.2)) +
    geom_line(position = position_dodge(0.2)) +
    scale_colour_manual(values = c(Female = "hotpink", Male = "blue")) +
    labs(title = metab,
         x = "Lifestyle Group", y = "Estimated Mean (±95% CI)", colour = "Sex") +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 9, face = "bold"), legend.position = "none",
      aspect.ratio = 1, panel.border = element_rect(color = "black", fill = NA))
})
# Display plots in 4×4 grid pages
for (i in seq(1, length(interaction_plots), by = 16)) {
  page <- interaction_plots[i:min(i + 15, length(interaction_plots))]
  print(ggarrange(plotlist = page, ncol = 4, nrow = 4))}


#####################################################################################################
# Checking for Lifestyle × Age interaction
#####################################################################################################
table(Data_for_NMR_ANALYSIS$lifestyle_group_MW)

# --- Evaluate models and store results (Lifestyle × Age) ---
lrt_summary_age <- map_dfr(metabolite_cols, function(metab) {
  tryCatch({
    df <- Data_for_NMR_ANALYSIS %>%
      dplyr::select(Metabolite = !!sym(metab), lifestyle_group_MW, Sex, Age) %>%
      drop_na()
    
    # Check factor levels / variation
    if (n_distinct(df$lifestyle_group_MW) < 2 || var(df$Age, na.rm = TRUE) == 0) {
      stop("Not enough information in lifestyle or age for contrasts")
    }
    
    # Reduced: no interaction
    reduced_model <- lm(Metabolite ~ lifestyle_group_MW + Sex + Age, data = df)
    # Full: with Lifestyle × Age interaction
    full_model    <- lm(Metabolite ~ lifestyle_group_MW * Age + Sex, data = df)
    
    lrt   <- anova(reduced_model, full_model)
    p_val <- lrt$`Pr(>F)`[2]
    
    tibble(
      Metabolite   = metab,
      LRT_p_age    = round(p_val, 6),
      Better_Model = ifelse(p_val < 0.05,
                            "Full (with L×Age interaction)",
                            "Reduced (no L×Age interaction)")
    )
  }, error = function(e) {
    message(glue::glue("Skipping {metab} (Age×Lifestyle) due to error: {e$message}"))
    NULL
  })
})

# --- Count how many times each model was better (Age × Lifestyle) ---
full_better_age    <- sum(lrt_summary_age$Better_Model == "Full (with L×Age interaction)")
reduced_better_age <- sum(lrt_summary_age$Better_Model == "Reduced (no L×Age interaction)")
total_age          <- full_better_age + reduced_better_age

cat("\nLifestyle × Age interaction summary:\n")
cat(paste("🧪 Full model (with L×Age) better in", full_better_age, "metabolites (",
          round(100 * full_better_age / total_age, 1), "% )\n"))
cat(paste("⚪ Reduced model (no L×Age) better in", reduced_better_age, "metabolites (",
          round(100 * reduced_better_age / total_age, 1), "% )\n"))

# --- Metabolites with evidence for Lifestyle × Age interaction ---
int_metabs_age <- lrt_summary_age %>%
  filter(Better_Model == "Full (with L×Age interaction)") %>%
  pull(Metabolite)

# --- Generate plots for significant Lifestyle × Age interactions ---
interaction_plots_age <- map(int_metabs_age, function(metab) {
  df <- Data_for_NMR_ANALYSIS %>%
    dplyr::select(Metabolite = !!sym(metab), lifestyle_group_MW, Sex, Age) %>%
    drop_na()
  
  # Fit model with Lifestyle × Age interaction
  fit <- lm(Metabolite ~ lifestyle_group_MW * Age + Sex, data = df)
  
  # Age range (5th–95th percentile) for smooth curves
  age_seq <- seq(
    quantile(df$Age, 0.05, na.rm = TRUE),
    quantile(df$Age, 0.95, na.rm = TRUE),
    length.out = 100
  )
  
  # Fix Sex at reference (e.g., Female) to visualize Lifestyle × Age
  ref_sex <- levels(df$Sex)[1]
  
  newd <- expand.grid(
    Age               = age_seq,
    lifestyle_group_MW = levels(df$lifestyle_group_MW),
    Sex               = factor(ref_sex, levels = levels(df$Sex))
  )
  
  pred <- as.data.frame(predict(fit, newdata = newd, interval = "confidence"))
  newd <- cbind(newd, pred)
  
  ggplot(newd,
         aes(x = Age, y = fit,
             color = lifestyle_group_MW,
             fill  = lifestyle_group_MW,
             group = lifestyle_group_MW)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr),
                alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    scale_color_manual(
      values = c(
        "Non_Urban" = "green",
        "Urban"     = "red"
      )
    ) +
    scale_fill_manual(
      values = c(
        "Non_Urban" = "green",
        "Urban"     = "red"
      )
    ) +
    labs(
      title    = metab,
      subtitle = paste("Lifestyle × Age interaction (Sex fixed at", ref_sex, ")"),
      x        = "Age",
      y        = "Estimated mean (±95% CI)"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title   = element_text(hjust = 0.5, size = 9, face = "bold"),
      legend.position = "none",
      aspect.ratio = 1,
      panel.border = element_rect(color = "black", fill = NA)
    )
})

# --- Display Lifestyle × Age interaction plots in 4×4 grid pages ---
for (i in seq(1, length(interaction_plots_age), by = 16)) {
  page <- interaction_plots_age[i:min(i + 15, length(interaction_plots_age))]
  print(ggarrange(plotlist = page, ncol = 4, nrow = 4))
}

