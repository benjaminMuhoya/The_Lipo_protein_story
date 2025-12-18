setwd("/Users/bm0211/RegEx/Water_Energy/2025_Metabolomic_Analysis/NMR_metabolomics/")
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
library(lme4)
library(broom.mixed)
library(lmerTest)
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
library(stringi)
# Load the data
## The combined file with those that pass QC - Whole Dataset
All_data <- import("/Users/bm0211/RegEx/Water_Energy/Scott_CO/THGP_database_TurkanaOnly_corrected_merged_2025-10-02.csv")
##Marina's APO-E list
APOE_Turks <- import("/Users/bm0211/RegEx/Water_Energy/2025_Metabolomic_Analysis/turkana_apoe_11Aug2025.txt")
head(APOE_Turks)
## The original NMR sample Map 
more_data_for_unmatched <- import("/Users/bm0211/RegEx/Water_Energy/HILIC_ms_analysis/NMR_Final_selection_sample_map.csv")
##Name of file from CVD_Risk Calculation Output ...Also contains Location data and coordinates
Cleaned_metadata_file <- import("/Users/bm0211/RegEx/Turkana_CVD/cvd_risk_results_male_female.csv") 
colnames(Cleaned_metadata_file)

Key_with_age_sex <- All_data %>%
  dplyr::select(Unique.ID, Nightingale_NMR_Sample.id, Age, Sex) %>%
  dplyr::filter(!is.na(Nightingale_NMR_Sample.id) & Nightingale_NMR_Sample.id != "")
tail(Key_with_age_sex)
dim(Key_with_age_sex)
head(more_data_for_unmatched)
# Add a new column indicating which UIDs are present in Key_with_age_sex
merged_df <- more_data_for_unmatched %>%
  mutate(Passed_QC = ifelse(UID %in% Key_with_age_sex$Unique.ID, 1, NA)) %>%
  dplyr::select(UID, Passed_QC)
head(merged_df)
merged_df$passedProcessing <- 1

# Check original
All_data <- All_data %>%
  mutate(Sex = ifelse(Sex %in% c("Male", "Female"), Sex, NA))
# Clean and create new columns
All_data <- All_data %>%
  mutate(
    medication_raw = If.you.are.on.any.medications.please.list.them.here,
    medication_clean = tolower(str_trim(medication_raw)),
    # Rename all hypertension-related mentions to "hypertension"
    medication_clean = case_when(
      grepl("htn medication|\\bbp\\b|\\bhtn\\b|hypertension", medication_clean, ignore.case = TRUE) ~ "hypertension",
      medication_clean == "" | is.na(medication_clean) ~ "no_medication",
      TRUE ~ medication_clean),
    # Create yes/no column preserving hypertension
    medication_yes_no = case_when(
      medication_clean == "hypertension" ~ "Hypertension",
      medication_clean == "no_medication" ~ "no",
      TRUE ~ "yes"))
# Check results
# Columns needed to calculate CVD Risk and composite score.
table(All_data$medication_clean)
table(All_data$Body.fat.percentage)
table(All_data$BMI)
table(All_data$Waist.circumference.1.cm.)
table(All_data$medication_yes_no)
table(All_data$Sex)
table(All_data$Age)
table(All_data$Total.cholesterol.mg.dL.)
table(All_data$HDL.cholesterol.mg.dL)
table(All_data$LDL.cholesterol.mg.dL)
table(All_data$Triglycerides.mg.dL)
table(All_data$Blood.pressure.mm.Hg.systolic.)
table(All_data$Blood.pressure.mm.Hg.diastolic.)
table(All_data$Tobacco.yes.no)
table(All_data$Type.I.diabetes)
table(All_data$Type.II.diabetes)
table(All_data$medication_yes_no)
table(All_data$Glucose.level.mg.dL.)
table(All_data$MW_scores_lifestyle)
write.csv(All_data, "/Users/bm0211/RegEx/Water_Energy/2025_Metabolomic_Analysis/RAW_Data_for_NMR_ANALYSIS.csv", row.names = FALSE)
getwd()
## run Framingham_CVD_turkana_script.py and import(RAW_Data_for_NMR_ANALYSIS_python_ran.csv)
RAW_Data_for_NMR_ANALYSIS <- import("/Users/bm0211/RegEx/Water_Energy/2025_Metabolomic_Analysis/RAW_Data_for_NMR_ANALYSIS_python_ran.csv")
colnames(RAW_Data_for_NMR_ANALYSIS)
dim(RAW_Data_for_NMR_ANALYSIS)
table(RAW_Data_for_NMR_ANALYSIS$MW_scores_lifestyle, useNA = "always")
RAW_Data_for_NMR_ANALYSIS <- RAW_Data_for_NMR_ANALYSIS %>%
  filter(!is.na(MW_scores_lifestyle), MW_scores_lifestyle != "")
# Ensure lifestyle_group is properly classified
RAW_Data_for_NMR_ANALYSIS <- RAW_Data_for_NMR_ANALYSIS %>%
  mutate(lifestyle_group_MW = ifelse(MW_scores_lifestyle %in% c("Pastoralist", "PeriUrban"),
                                     "Non_Urban", "Urban"))
table(RAW_Data_for_NMR_ANALYSIS$lifestyle_group_MW)
##
ggplot(RAW_Data_for_NMR_ANALYSIS, aes(x = Framingham_10yr_CVD_Risk)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
  labs(title = "Framingham 10-Year CVD Risk", x = "CVD Risk (%)", y = "Count") + theme_minimal()
###############################################################################
###############################################################################
###############################################################################
###############################################################################
# Define clinical thresholds
thresholds <- list(
  BMI = 28,
  Waist_circumference = 100,                 # You can adjust based on sex if needed
  Body_fat = 25,                             # Approximate general threshold (adjust if sex-specific needed)
  Total_cholesterol = 200,
  HDL = 40,
  LDL = 130,
  Triglycerides = 150,
  Glucose = 100,
  Systolic_BP = 130,
  Diastolic_BP = 80)
###############################################################################
###############################################################################
# Compute composite score
RAW_Data_for_NMR_ANALYSIS$composite_score <- apply(RAW_Data_for_NMR_ANALYSIS, 1, function(row) {
  biomarker_flags <- c(
    as.numeric(row["BMI"]) > thresholds$BMI,
    as.numeric(row["Waist.circumference.1.cm."]) > thresholds$Waist_circumference,
    as.numeric(row["Body.fat.percentage"]) > thresholds$Body_fat,
    as.numeric(row["Total.cholesterol.mg.dL."]) > thresholds$Total_cholesterol,
    as.numeric(row["HDL.cholesterol.mg.dL"]) < thresholds$HDL,
    as.numeric(row["LDL.cholesterol.mg.dL"]) > thresholds$LDL,
    as.numeric(row["Triglycerides.mg.dL"]) > thresholds$Triglycerides,
    as.numeric(row["Glucose.level.mg.dL."]) > thresholds$Glucose,
    as.numeric(row["Blood.pressure.mm.Hg.systolic."]) > thresholds$Systolic_BP,
    as.numeric(row["Blood.pressure.mm.Hg.diastolic."]) > thresholds$Diastolic_BP
  )
  # Count valid (non-NA) biomarkers
  num_non_missing <- sum(!is.na(biomarker_flags))
  if (num_non_missing >= 3) {
    return(sum(biomarker_flags, na.rm = TRUE) / num_non_missing)
  } else {
    return(NA)
  }
})
table(RAW_Data_for_NMR_ANALYSIS$composite_score)
##
ggplot(RAW_Data_for_NMR_ANALYSIS, aes(x = composite_score)) +
  geom_histogram(binwidth = 0.05, fill = "blue4", color = "black", boundary = 0) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
  labs(title = "Distribution of Composite Cardiometabolic Risk Score",
       x = "Proportion of Biomarkers Above Clinical Cutoff", y = "Number of Individuals") + theme_minimal()
###############################################################################
RAW_Data_for_NMR_ANALYSIS <- RAW_Data_for_NMR_ANALYSIS %>%
  dplyr::filter(!is.na(MW_scores_lifestyle) & MW_scores_lifestyle != "")
# Calculate chi-square
table(RAW_Data_for_NMR_ANALYSIS$MW_scores_lifestyle)
dim(RAW_Data_for_NMR_ANALYSIS)
table_summary_composite <- table(RAW_Data_for_NMR_ANALYSIS$MW_scores_lifestyle,
                                 RAW_Data_for_NMR_ANALYSIS$composite_score)
table(RAW_Data_for_NMR_ANALYSIS$lifestyle_group_MW)
###############################################################################
model_lifestyle <- glm(composite_score ~ MW_scores_lifestyle + Age + Sex,
                       family = quasibinomial(link = "logit"),
                       data = RAW_Data_for_NMR_ANALYSIS)
summary(model_lifestyle)
# Estimated marginal means on response scale
em_lifestyle <- emmeans(model_lifestyle, pairwise ~ MW_scores_lifestyle, type = "response")
summary(em_lifestyle)
# Extract emmeans estimates
emm_df <- as.data.frame(em_lifestyle$emmeans)
# Extract odds ratios for each pairwise comparison
contrast_df <- as.data.frame(em_lifestyle$contrasts)
# Custom labels for odds ratios
contrast_labels <- data.frame(
  x_start = c(1, 1, 2), x_end   = c(2, 3, 3),
  y       = c(0.20, 0.22, 0.24),  # adjust based on your ymax
  label   = c("OR = 1.03, p = 0.94",
              "OR = 0.67, p < 0.001",
              "OR = 0.65, p < 0.001"))

# Define the color mapping
lifestyle_colors <- c("PeriUrban" = "green", "Urban" = "#FC8D62", "Pastoralist" = "#8DA0CB")
ggplot(emm_df, aes(x = MW_scores_lifestyle, y = prob, fill = MW_scores_lifestyle)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  geom_text(aes(label = sprintf("%.3f", prob)), vjust = -3, size = 4.5) +
  # Add odds ratio brackets
  geom_segment(data = contrast_labels,
               aes(x = x_start, xend = x_end, y = y, yend = y),
               inherit.aes = FALSE) +
  geom_text(data = contrast_labels,
            aes(x = (x_start + x_end) / 2, y = y + 0.005, label = label),
            size = 4, vjust = 0, inherit.aes = FALSE) +
  scale_fill_manual(values = lifestyle_colors) +
  labs(title = "CVD Risk by Lifestyle - Lea 2020 COMPOSITE SCORE",
       x = "Lifestyle Group", y = "Proportion of Risky Biomarkers",
       fill = "Lifestyle") + theme_minimal(base_size = 14) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        legend.position = "none"  # remove legend if not needed
  ) + ylim(0, max(contrast_labels$y) + 0.05)

###############################################################################
###############################################################################
###############################################################################
ggplot(RAW_Data_for_NMR_ANALYSIS, aes(x = composite_score, fill = MW_scores_lifestyle)) +
  geom_histogram(binwidth = 0.05, color = "black", boundary = 0) +
  facet_wrap(~ MW_scores_lifestyle, nrow = 1) +
  scale_fill_manual(values = lifestyle_colors) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
  labs( title = "Composite Score Distribution by Lifestyle Group",
        x = "Proportion of Biomarkers Above Clinical Cutoff",
        y = "Number of Individuals", fill = "Lifestyle" ) + theme_minimal()
########
ggplot(RAW_Data_for_NMR_ANALYSIS, aes(x = composite_score, fill = MW_scores_lifestyle)) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values = lifestyle_colors) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
  labs(title = "Density of COMPOSITE Score by Lifestyle Group",
       x = "Proportion of Biomarkers Above Clinical Cutoff",
       y = "Density", fill = "Lifestyle", color = "Lifestyle") + theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        legend.position = "none"  # remove legend if not needed
  ) 
###############################################################################
RAW_Data_for_NMR_ANALYSIS <- RAW_Data_for_NMR_ANALYSIS %>%
  dplyr::filter(!is.na(MW_scores_lifestyle) & MW_scores_lifestyle != "")
###############################################################################
# Framingham
model_framingham <- lm(Framingham_10yr_CVD_Risk ~ MW_scores_lifestyle + Age + Sex,
                       data = RAW_Data_for_NMR_ANALYSIS)
summary(model_framingham)
########
em_framingham <- emmeans(model_framingham, pairwise ~ MW_scores_lifestyle)
emm_df_fram <- as.data.frame(em_framingham$emmeans)
contrast_df_fram <- as.data.frame(em_framingham$contrasts)
######
# Parse group names from the contrast column (e.g., "A - B")
group_names <- strsplit(as.character(contrast_df_fram$contrast), " - ")
x_start_vals <- sapply(group_names, function(g) match(g[1], levels(emm_df_fram$MW_scores_lifestyle)))
x_end_vals   <- sapply(group_names, function(g) match(g[2], levels(emm_df_fram$MW_scores_lifestyle)))
# Format the p-value dynamically
format_p <- function(p) {
  if (p < 0.001) {
    return("p < 0.001")
  } else {
    return(paste0("p = ", signif(p, 2)))
  }
}
# Generate clean labels
contrast_labels_fram <- data.frame(
  x_start = x_start_vals, x_end   = x_end_vals,
  y       = seq(from = max(emm_df_fram$emmean) + 1.5, by = 0.7, length.out = nrow(contrast_df_fram)),
  label   = paste0( gsub(" - ", " − ", contrast_df_fram$contrast), ": ",
                    ifelse(contrast_df_fram$estimate > 0, "+", ""),
                    round(contrast_df_fram$estimate, 2), ", ",
                    sapply(contrast_df_fram$p.value, format_p)))
#######
ggplot(emm_df_fram, aes(x = MW_scores_lifestyle, y = emmean, fill = MW_scores_lifestyle)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  geom_text(aes(label = sprintf("%.1f", emmean)), vjust = -3, size = 4.5) +
  geom_segment(data = contrast_labels_fram,
               aes(x = x_start, xend = x_end, y = y, yend = y),
               inherit.aes = FALSE) +
  geom_text(data = contrast_labels_fram,
            aes(x = (x_start + x_end) / 2, y = y + 0.4, label = label),
            size = 4, vjust = 0, inherit.aes = FALSE) +
  scale_fill_manual(values = lifestyle_colors) +
  labs(title = "Estimated Framingham 10-Year CVD Risk by Lifestyle Group",
       x = "Lifestyle Group",y = "Predicted Risk Score (%)",
       fill = "Lifestyle") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        legend.position = "none"
  ) + ylim(0, max(contrast_labels_fram$y) + 2)
# Density plot of Framingham CVD Risk
RAW_Data_for_NMR_ANALYSIS %>%
  filter(!is.na(Framingham_10yr_CVD_Risk), !is.na(MW_scores_lifestyle)) %>%
  ggplot(aes(x = Framingham_10yr_CVD_Risk, fill = MW_scores_lifestyle)) +
  geom_density(alpha = 0.3, adjust = 1.2) +
  scale_fill_manual(values = lifestyle_colors) +
  scale_x_continuous(limits = c(0, 40), breaks = seq(0, 40, 5)) +
  labs(title = "Framingham 10-Year CVD Risk by Lifestyle Group",
       x = "CVD Risk Score (%)", y = "Density", fill = "Lifestyle") +
  theme_minimal() + theme(panel.grid.major.x = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.background = element_rect(fill = "white", color = NA),
                          plot.background = element_rect(fill = "white", color = NA),
                          legend.position = "none")