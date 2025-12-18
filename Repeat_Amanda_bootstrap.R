setwd("/Users/bm0211/RegEx/Water_Energy/2025_Metabolomic_Analysis/NMR_metabolomics/")
# Load necessary libraries
library(effsize)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(rio)
library(gridExtra)
library(grid) 
library(reshape2)
library(tidyverse)
library(tidyr)
library(broom)
library(broom.mixed)
library(ggpubr)
library(ggfortify)
library(patchwork)
library(wesanderson)
library(emmeans)
library(ggvenn)
library(UpSetR)
library(dplyr)
library(readr)
library(stringr)
library(VennDiagram)
library(scales)
library(mediation)
library(leaflet)
library(leaflet.minicharts)
library(htmlwidgets)
library(fuzzyjoin)
# Load the data
# This data is output by the previous script 
# Framingham_CVD_Risk_and_CLEAN_Data_for_downstream.R
RAW_Data_for_NMR_ANALYSIS <- import("/Users/bm0211/RegEx/Water_Energy/2025_Metabolomic_Analysis/RAW_Data_for_NMR_ANALYSIS.csv")
colnames(RAW_Data_for_NMR_ANALYSIS)
##Subset data to only the needed columns for further NMR Metabolomics
columns_to_keep <- c(1,25,31,590, 49,600,602,875, 888,890,891, 892, 903, 905, 625:874, 589,489:509)
# Subset the data to include only the selected columns
RAW_Data_for_NMR_ANALYSIS <- RAW_Data_for_NMR_ANALYSIS %>% dplyr::select(all_of(columns_to_keep))
dim(RAW_Data_for_NMR_ANALYSIS)
colnames(RAW_Data_for_NMR_ANALYSIS)
table(RAW_Data_for_NMR_ANALYSIS$MW_scores_lifestyle)
table(RAW_Data_for_NMR_ANALYSIS$Sex)
# Rename columns by removing "Nightingale_NMR_" prefix
colnames(RAW_Data_for_NMR_ANALYSIS) <- gsub("^Nightingale_NMR_", "", colnames(RAW_Data_for_NMR_ANALYSIS))
# Remove those ending with _pct & .PCT i.e. Percentages. 
RAW_Data_for_NMR_ANALYSIS <- RAW_Data_for_NMR_ANALYSIS %>%
  dplyr::select(-matches("_pct$"),  -matches("(?i)\\.PCT$"))
#####################################################################################################
#####################################################################################################
# -------------------------------- Process the columns for analysis ---------------------------------
#####################################################################################################
#####################################################################################################
colnames(RAW_Data_for_NMR_ANALYSIS)
# Process all the Numeric ones - scale and remove outliers
columns_to_process <- c(2,6:11, 13:64,67:192, 195:204) # Age, TL_108 , Diet Metrics included ...I excluded Ratios and factor
# Factor Columns
RAW_Data_for_NMR_ANALYSIS$Sex <- as.factor(RAW_Data_for_NMR_ANALYSIS$Sex)
table(RAW_Data_for_NMR_ANALYSIS$Sex)
RAW_Data_for_NMR_ANALYSIS <- RAW_Data_for_NMR_ANALYSIS %>%
  mutate(lifestyle_group_MW = recode(MW_scores_lifestyle,
                                "Pastoralist" = "Non_Urban",
                                "PeriUrban"   = "Non_Urban",
                                "Urban"       = "Urban") |> 
      factor(levels = c("Non_Urban", "Urban")))
RAW_Data_for_NMR_ANALYSIS$lifestyle_group_MW <- as.factor(RAW_Data_for_NMR_ANALYSIS$lifestyle_group_MW)
table(RAW_Data_for_NMR_ANALYSIS$lifestyle_group_MW)
#####################################################################################################
# Function to convert columns to numeric
convert_to_numeric <- function(df, columns) {
  for (col in columns) {
    column_name <- colnames(df)[col]
    df[[column_name]] <- as.numeric(as.character(df[[column_name]]))}
  return(df)}
RAW_Data_for_NMR_ANALYSIS <- convert_to_numeric(RAW_Data_for_NMR_ANALYSIS, columns_to_process)
# -------- Function to remove outliers RAW data -------------
remove_outliers_iqr <- function(df, columns) {
  for (col in columns) {
    column_name <- colnames(df)[col]
    if (is.numeric(df[[column_name]])) {
      Q1 <- quantile(df[[column_name]], 0.25, na.rm = TRUE)
      Q3 <- quantile(df[[column_name]], 0.75, na.rm = TRUE)
      IQR <- Q3 - Q1
      lower_bound <- Q1 - 1.5 * IQR
      upper_bound <- Q3 + 1.5 * IQR
      # Replace outliers with NA
      df[[column_name]][df[[column_name]] < lower_bound | df[[column_name]] > upper_bound] <- NA
    }
  }
  return(df)
}
RAW_Data_for_NMR_ANALYSIS <- remove_outliers_iqr(RAW_Data_for_NMR_ANALYSIS, columns_to_process)

# Standardize and normalize columns
standardize_columns <- function(df, columns) {
  for (col in columns) {
    column_name <- colnames(df)[col]
    log_values <- log10(df[[column_name]])
    df[[column_name]] <- scale(log_values, center = TRUE, scale = TRUE)
  }
  return(df)
}
Data_for_NMR_ANALYSIS <- standardize_columns(RAW_Data_for_NMR_ANALYSIS, columns_to_process)
table(Data_for_NMR_ANALYSIS$lifestyle_group_MW)
# -------- Function to remove outliers and replace with NA (transformed data) -------------
remove_outliers_iqr <- function(df, columns) {
  for (col in columns) {
    column_name <- colnames(df)[col]
    if (is.numeric(df[[column_name]])) {
      Q1 <- quantile(df[[column_name]], 0.25, na.rm = TRUE)
      Q3 <- quantile(df[[column_name]], 0.75, na.rm = TRUE)
      IQR <- Q3 - Q1
      lower_bound <- Q1 - 1.5 * IQR
      upper_bound <- Q3 + 1.5 * IQR
      # Replace outliers with NA
      df[[column_name]][df[[column_name]] < lower_bound | df[[column_name]] > upper_bound] <- NA
    }
  }
  return(df)
}
Data_for_NMR_ANALYSIS <- remove_outliers_iqr(Data_for_NMR_ANALYSIS, columns_to_process)
colnames(Data_for_NMR_ANALYSIS)
#####################################################################################################
# Checking for Normality
# Define metabolite columns
metabolite_cols <- names(Data_for_NMR_ANALYSIS[, c(9:204)]) # everything including inflammation ones
# Drop NA from batch column
Data_for_NMR_ANALYSIS_filtered <- Data_for_NMR_ANALYSIS
# ---- Run Shapiro test and store results ----
normality_results <- data.frame()
for (metab in metabolite_cols) {
  current_data <- Data_for_NMR_ANALYSIS_filtered %>%
    dplyr::select(all_of(metab)) %>%
    drop_na()
  values <- current_data[[1]]
  if (length(values) >= 4 && length(unique(values)) > 1) {
    shapiro <- shapiro.test(values)
    
    normality_results <- rbind(normality_results, data.frame(
      Metabolite = metab,
      N = length(values),
      W_statistic = as.numeric(shapiro$statistic),
      p_value = as.numeric(shapiro$p.value),
      Deviates_from_normal = ifelse(shapiro$p.value < 0.05, "Yes", "No")
    ))
  } else {
    reason <- if (length(values) < 4) "Insufficient data" else "Identical values"
    normality_results <- rbind(normality_results, data.frame(
      Metabolite = metab,
      N = length(values),
      W_statistic = NA,
      p_value = NA,
      Deviates_from_normal = reason
    ))
  }
}
# ---- Filter to usable results and order by deviation ----
plot_data <- normality_results %>%
  filter(Deviates_from_normal %in% c("Yes", "No")) %>%
  mutate(order_priority = ifelse(Deviates_from_normal == "Yes", 1, 2)) %>%
  arrange(order_priority, p_value)
# ---- Create annotated histograms ----
pdf("Metabolite_Histograms_Annotated_LOG.pdf", width = 14, height = 10)
plot_list <- list()
count <- 0
for (i in seq_len(nrow(plot_data))) {
  metab <- plot_data$Metabolite[i]
  values <- Data_for_NMR_ANALYSIS_filtered[[metab]]
  values <- values[!is.na(values)]
  # Build label
  label_text <- paste0(
    "N = ", plot_data$N[i],
    ", W = ", signif(plot_data$W_statistic[i], 4),
    ", p = ", signif(plot_data$p_value[i], 3),
    "\nDeviates from Normal (Shapiro.test)? ", plot_data$Deviates_from_normal[i])
  
  p <- ggplot(data.frame(x = values), aes(x = x)) +
    geom_histogram(bins = 30, fill = "skyblue", color = "black") +
    labs(title = metab, subtitle = label_text, x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(plot.subtitle = element_text(size = 8))
  plot_list[[length(plot_list) + 1]] <- p
  count <- count + 1
  if (count %% 25 == 0) {
    do.call("grid.arrange", c(plot_list, ncol = 5))
    plot_list <- list()
  }
}
# Print final page
if (length(plot_list) > 0) {
  do.call("grid.arrange", c(plot_list, ncol = 5))
}
dev.off()

#####################################################################################################
# Define metabolite columns
metabolite_cols <- names(Data_for_NMR_ANALYSIS[, c(10,11, 15:86, 89:198)]) ##Exclude inflammation columns
inflammation_columns <-names(Data_for_NMR_ANALYSIS[, c(9,12:14,87,88,188,199:204)])
diet_columns <- names(Data_for_NMR_ANALYSIS[, c(6,7,57:66)])
# Define non-normal metabolites to skip
non_normal_metabolites <- c(
  "XXL_VLDL_P", "XXL_VLDL_L", "XXL_VLDL_PL", "XXL_VLDL_C", "XXL_VLDL_CE",
  "XXL_VLDL_FC", "XXL_VLDL_TG", "XL_VLDL_P", "XL_VLDL_L", "XL_VLDL_PL",
  "XL_VLDL_FC", "XL_VLDL_TG")
# Filter to only normal metabolites
normal_metabolites <- setdiff(metabolite_cols, non_normal_metabolites)
## Factor Columns
Data_for_NMR_ANALYSIS$Sex <- as.factor(Data_for_NMR_ANALYSIS$Sex)
# lock factor order globally
Data_for_NMR_ANALYSIS$lifestyle_group_MW <- factor(
  Data_for_NMR_ANALYSIS$lifestyle_group_MW,
  levels = c("Non_Urban", "Urban")
)
table(Data_for_NMR_ANALYSIS$lifestyle_group_MW)
Data_for_NMR_ANALYSIS$TA_score_TL_108 <- as.numeric(Data_for_NMR_ANALYSIS$TA_score_TL_108)
#####################################################################################################
## Repeat Amanda's 2020 Analysis
####################################################################################################
## CLEAN + SIMPLIFIED BOOTSTRAP PIPELINE FOR LIFESTYLE / SEX / AGE EFFECTS
## Purpose:
##   1. Compute bootstrap estimates (β ± CI, p<0.05 proportion) for LDL, HDL, Chol and Trig
##   2. Display one combined bar plot showing effects of each predictor on each biomarker
##   3. Generate interaction plots only for significant interaction terms
## Notes:
##  - Sampling uses balanced N for each group ..500N Urban and 500N Non_Urban 
####################################################################################################
set.seed(42)
## ── Utility functions ───────────────────────────────────────────────────────────
# Effect size (Cohen’s d) for two-group comparisons
cohens_d_2group <- function(y, g) {
  g <- droplevels(g)
  x <- y[g == levels(g)[1]]; z <- y[g == levels(g)[2]]
  sp <- sqrt(((length(x)-1)*var(x) + (length(z)-1)*var(z)) / (length(x)+length(z)-2))
  if (!is.finite(sp) || sp == 0) return(NA_real_)
  (mean(z) - mean(x)) / sp
}

# Keep non-outliers using IQR rule
iqr_keep <- function(v) {
  q1 <- quantile(v, 0.25, na.rm = TRUE); q3 <- quantile(v, 0.75, na.rm = TRUE)
  i  <- IQR(v, na.rm = TRUE)
  (v >= q1 - 1.5*i) & (v <= q3 + 1.5*i)
}


## ── Bootstrap for any linear model formula ──────────────────────────────────────
bootstrap_formula_one <- function(bm, dat, form, coef_keys,
                                  lifestyle_levels = c("Non_Urban","Urban"),
                                  sex_levels = c("Female","Male"),
                                  B = 10000, n_target = 1000) {
  # Build filtered dataset
  base <- dat %>%
    transmute(
      y = as.numeric(.data[[bm]]),
      Lifestyle = factor(.data[["lifestyle_group_MW"]], levels = lifestyle_levels),
      Sex       = factor(.data[["Sex"]],               levels = sex_levels),
      Age       = as.numeric(.data[["Age"]])
    ) %>%
    filter(is.finite(y), is.finite(Age)) %>%
    filter(iqr_keep(y))
  
  # Skip biomarkers with insufficient data
  if (nrow(base) < 50 || var(base$y) == 0 ||
      nlevels(base$Lifestyle) < 2 || nlevels(base$Sex) < 2) {
    return(tibble(
      Biomarker = bm, N_available = nrow(base),
      Model = deparse(form),
      Term = names(coef_keys),
      beta_median = NA_real_, beta_L = NA_real_, beta_U = NA_real_,
      prop_p_lt_0.05 = NA_real_
    ))
  }
  # ---- Resample proportional to the original strata mix (Lifestyle × Sex) ----
  base$Stratum <- interaction(base$Lifestyle, base$Sex, drop = TRUE)
  tab <- table(base$Stratum)               # observed counts per stratum
  str_names <- names(tab)
  
  # Allocate n_target across strata in proportion to observed frequencies
  n_per <- as.integer(round(tab / sum(tab) * n_target))
  
  # Ensure at least 1 from any existing stratum (for model stability)
  n_per <- pmax(1L, n_per)
  
  # Fix rounding so totals sum exactly to n_target
  diff_total <- n_target - sum(n_per)
  if (diff_total != 0) {
    j <- which.max(tab)                    # adjust the largest stratum
    n_per[j] <- n_per[j] + diff_total
  }
  names(n_per) <- str_names
  
  # Bootstrap loop
  out <- replicate(B, {
    idx <- unlist(lapply(names(n_per), function(slv) {
      pool <- which(base$Stratum == slv); size <- n_per[slv]
      sample(pool, size = size, replace = size > length(pool))
    }))
    s <- base[idx, , drop = FALSE]
    if (var(s$y) == 0 ||
        nlevels(droplevels(s$Lifestyle)) < 2 ||
        nlevels(droplevels(s$Sex)) < 2) {
      v <- rep(NA_real_, length(coef_keys)*2)
      names(v) <- c(names(coef_keys), paste0("p_", names(coef_keys)))
      return(v)
    }
    
    fit <- lm(form, data = s)
    co <- summary(fit)$coefficients
    
    betas <- sapply(coef_keys, function(nm) if (nm %in% rownames(co)) co[nm, "Estimate"] else NA_real_)
    pvals <- sapply(coef_keys, function(nm) if (nm %in% rownames(co)) co[nm, "Pr(>|t|)"] else NA_real_)
    names(betas) <- names(coef_keys)
    names(pvals) <- paste0("p_", names(coef_keys))
    c(betas, pvals)
  })
  
  bx <- out[names(coef_keys), , drop = FALSE]
  px <- out[paste0("p_", names(coef_keys)), , drop = FALSE]
  
  tibble(
    Biomarker   = bm,
    N_available = nrow(base),
    Model       = deparse(form),
    Term        = names(coef_keys),
    beta_median = apply(bx, 1, function(v) median(v, na.rm = TRUE)),
    beta_L      = apply(bx, 1, function(v) quantile(v, 0.025, na.rm = TRUE)),
    beta_U      = apply(bx, 1, function(v) quantile(v, 0.975, na.rm = TRUE)),
    prop_p_lt_0.05 = sapply(rownames(px), function(rn) mean(px[rn, ] < 0.05, na.rm = TRUE))
  )
}

# Wrapper to run bootstrap for all biomarkers
bootstrap_formula_all <- function(biomarkers, dat, form, coef_keys,
                                  B = 10000, n_target = 1000) {
  do.call(rbind, lapply(
    biomarkers,
    bootstrap_formula_one,
    dat = dat, form = form, coef_keys = coef_keys,
    B = B, n_target = n_target
  ))
}

## ── Define models ───────────────────────────────────────────────────────────────
form_base <- y ~ Lifestyle + Sex + Age
coef_base <- c("Lifestyle (Urban−Non_Urban)" = "LifestyleUrban",
               "Sex (Male−Female)" = "SexMale",
               "Age (per unit)" = "Age")

form_LxA <- y ~ Lifestyle * Age + Sex
coef_LxA <- c(coef_base, "Lifestyle×Age" = "LifestyleUrban:Age")

form_LxS <- y ~ Lifestyle * Sex + Age
coef_LxS <- c(coef_base, "Lifestyle×Sex" = "LifestyleUrban:SexMale")

forms <- list(
  Base = list(form=form_base, coef=coef_base),
  LxA  = list(form=form_LxA,  coef=coef_LxA),
  LxS  = list(form=form_LxS,  coef=coef_LxS)
)

## ── Run bootstrap for all models ────────────────────────────────────────────────
colnames(Data_for_NMR_ANALYSIS)
TO_BOOT_STRAP <- Data_for_NMR_ANALYSIS %>%
  mutate(lifestyle_group_MW = factor(lifestyle_group_MW, levels = c("Non_Urban","Urban")),
         Sex = factor(Sex, levels = c("Female","Male")),
         Age = as.numeric(Age))

biomarkers <- c("Total.cholesterol.mg.dL.","HDL.cholesterol.mg.dL",
                "LDL.cholesterol.mg.dL","Triglycerides.mg.dL",
                "COLO_APOLIPOPROTEINS_A1_g/L", "COLO_APOLIPOPROTEIN_B_g/L",
                "Non.HDL.mg.dL", "LBP_ELISA_Adjusted_dilutionFactor_finalConcentration_ug_ml",
                "S_HDL_P", "M_HDL_P", "L_HDL_P", "XL_HDL_P",
                "COLO_GLOBULIN_g/L", "COLO_C_REACTIVE_PROTEIN_mg/L",
                "COLO_ALBUMIN_SERUM_g/L")
#biomarkers <- names(TO_BOOT_STRAP)[9:204]
all_models <- bind_rows(lapply(names(forms), function(m){
  bootstrap_formula_all(biomarkers, TO_BOOT_STRAP,
                        forms[[m]]$form, forms[[m]]$coef,
                        B = 10000, n_target = 1000) %>%
    mutate(ModelTag = m)
}))

write_csv(all_models, "Bootstrap_AllModels_Lifestyle_Sex_Age_and_Interactions.csv")
head(all_models)

## ── Combined bar plot: predictors × biomarkers ──────────────────────────────────
plot_terms_all <- function(tab) {
  tab <- tab |>
    dplyr::distinct(Biomarker, Term, .keep_all = TRUE) |>
    dplyr::mutate(
      Sig  = ifelse(prop_p_lt_0.05 >= 0.95, "*", ""),
      Term = factor(Term, levels = c(
        "Lifestyle (Urban−Non_Urban)", "Sex (Male−Female)",
        "Age (per unit)", "Lifestyle×Age", "Lifestyle×Sex"
      ))
    )
  
  span <- with(tab, max(beta_U, na.rm=TRUE) - min(beta_L, na.rm=TRUE))
  
  ggplot(tab, aes(y = Term)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey55") +
    ggplot2::geom_errorbarh(aes(xmin = beta_L, xmax = beta_U), height = 0, linewidth = 1) +
    geom_point(aes(x = beta_median), size = 2.8) +
    # asterisk to the right of CI
    geom_text(aes(x = beta_U + 0.03*span, label = Sig), fontface = "bold", size = 5, na.rm = TRUE) +
    facet_wrap(~ Biomarker, ncol = 2, scales = "free_x") +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.12))) +
    labs(title = "Predictor effects across biomarkers (N=1000, 10,000 bootstraps)",
         subtitle = "Point = median β; line = 95% CI; asterisk = ≥95% of bootstraps p<0.05",
         x = "Median β (95% CI)", y = NULL) +
    theme_minimal(base_size = 14) +
    theme(strip.text = element_text(face = "bold", size = 13),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank()) + 
    coord_cartesian(clip = "off") +  # keeps asterisks visible past edges
    theme(panel.border = element_rect(color = "grey55", fill = NA, linewidth = 0.6),
          panel.spacing = unit(10, "pt"),
          strip.background = element_rect(fill = "grey96", color = "grey70")
    )
}

p_all <- plot_terms_all(all_models)
print(p_all)

# Save multi-page PDF with 9 biomarkers per page
outdir <- "Bootstrap_FIGURES"
dir.create(outdir, showWarnings = FALSE)

bio_chunks <- split(unique(all_models$Biomarker),
                    ceiling(seq_along(unique(all_models$Biomarker)) / 6))

pdf(file.path(outdir, "Predictor_effects_across_biomarkers.pdf"), width = 10, height = 8)
for (b in bio_chunks) {
  p_sub <- plot_terms_all(all_models %>% filter(Biomarker %in% b))
  print(p_sub)
}
dev.off()

## ── Identify significant interaction terms ─────────────────────────────────────
THR <- 0.95
sig_LxA <- all_models %>%
  filter(ModelTag == "LxA", Term == "Lifestyle×Age",
         is.finite(prop_p_lt_0.05), prop_p_lt_0.05 >= THR) %>%
  pull(Biomarker) %>% unique()

sig_LxS <- all_models %>%
  filter(ModelTag == "LxS", Term == "Lifestyle×Sex",
         is.finite(prop_p_lt_0.05), prop_p_lt_0.05 >= THR) %>%
  pull(Biomarker) %>% unique()

## ── clean df before plotting interactions ────────────────────────────
prep_bm_df <- function(bm, dat) {
  dat %>%
    transmute(
      y = as.numeric(.data[[bm]]),
      Lifestyle = factor(.data[["lifestyle_group_MW"]], levels = c("Non_Urban","Urban")),
      Sex = factor(.data[["Sex"]], levels = c("Female","Male")),
      Age = as.numeric(.data[["Age"]])
    ) %>%
    filter(is.finite(y), is.finite(Age)) %>%
    filter(iqr_keep(y))
}

## ── Interaction plots for significant biomarkers ────────────────────────────────
plot_LxA <- function(bm, dat, by_sex = "Female") {
  df <- prep_bm_df(bm, dat)
  if (nrow(df) < 50 || var(df$y) == 0) return(NULL)
  fit <- lm(y ~ Lifestyle * Age + Sex, data = df)
  
  a_seq <- seq(quantile(df$Age, 0.05, na.rm = TRUE),
               quantile(df$Age, 0.95, na.rm = TRUE), length.out = 100)
  newd <- expand.grid(Age = a_seq,
                      Lifestyle = levels(df$Lifestyle),
                      Sex = factor(by_sex, levels = levels(df$Sex)))
  pr <- as.data.frame(predict(fit, newd, interval = "confidence"))
  newd <- cbind(newd, pr)
  
  ggplot(newd, aes(x = Age, y = fit, color = Lifestyle, fill = Lifestyle)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    labs(title = paste0(bm, " — Lifestyle × Age interaction"),
         x = "Age", y = "Predicted trait") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "top")
}

plot_LxS <- function(bm, dat) {
  df <- prep_bm_df(bm, dat)
  if (nrow(df) < 50 || var(df$y) == 0) return(NULL)
  fit <- lm(y ~ Lifestyle * Sex + Age, data = df)
  A <- median(df$Age, na.rm = TRUE)
  newd <- expand.grid(Lifestyle = levels(df$Lifestyle),
                      Sex = levels(df$Sex),
                      Age = A)
  pr <- as.data.frame(predict(fit, newd, interval = "confidence"))
  newd <- cbind(newd, pr)
  
  ggplot(newd, aes(x = Sex, y = fit, color = Lifestyle, group = Lifestyle)) +
    geom_point(position = position_dodge(width = 0.4), size = 2.8) +
    geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.15,
                  position = position_dodge(width = 0.4)) +
    geom_line(position = position_dodge(width = 0.4)) +
    labs(title = paste0(bm, " — Lifestyle × Sex interaction"),
         subtitle = paste("Age fixed at", round(A,2)),
         x = "Sex", y = "Predicted trait") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "top")
}

plots_LxA <- map(sig_LxA, ~ plot_LxA(.x, TO_BOOT_STRAP))
plots_LxS <- map(sig_LxS, ~ plot_LxS(.x, TO_BOOT_STRAP))

# Save interactions to multi-page PDF (6 per page)
pdf(file.path(outdir, "Interaction_plots.pdf"), width = 8, height = 6)
plotsA <- plots_LxA[!sapply(plots_LxA, is.null)]
plotsS <- plots_LxS[!sapply(plots_LxS, is.null)]

for (i in seq_along(plotsA)) {
  print(plotsA[[i]])
  if (i %% 6 == 0 && i < length(plotsA)) grid::grid.newpage()
}

for (i in seq_along(plotsS)) {
  print(plotsS[[i]])
  if (i %% 6 == 0 && i < length(plotsS)) grid::grid.newpage()
}
dev.off()


####################################################################################################