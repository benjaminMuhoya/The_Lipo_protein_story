# Load necessary libraries
library(dplyr)
library(cluster)
library(factoextra)
library(readr)
library(tidyverse)
library(gridExtra) 
library(rio)
library(stringr)
library(ggplot2)
library(readxl)
library(stringr)
library(purrr)
library(patchwork)
library(gtools)
library(scales)
library(broom)
library(forcats)
# Set working directory (if needed)
setwd("/Users/Bm0211/RegEx/Water_Energy/")
###############################
####################
#########
# Define a vector with 20 distinct colors
my_colors_P <- c(
  "#FAD510", # Mustard Yellow
  "#008080", # Teal
  "lightpink", # Pink
  "green", # Green
  "#FF2400", # Scarlet Red
  "cyan", # Cyan
  "#90A959", # Olive Green
  "#9D858D", # Mauve
  "#A4243B", # Deep Carmine
  "#6495ED", # Cornflower Blue
  "#5B1A18", # Maroon
  "#FF00FF", # Magenta
  "#C6C8CA", # Light Gray
  "#E9A368", # Salmon Orange
  "#81A3A7", # Steel Blue
  "blue", # Gold
  "#D2691E", # Chocolate Brown
  "#4B0082", # Indigo
  "white", #
  "#7FFF00",  # Chartreuse
  "purple",
  "orange",
  "black", # lightgrey
  "darkgreen", 
  "#FF69B4", 
  
  "#00CED1", "#FFD700" 
)
# Step 1: Import datasets
HDL_Meta_analysis <- import("/Users/Bm0211/RegEx/Water_Energy/Cholesterol_Meta_analysis.csv")
colnames(HDL_Meta_analysis)
dim(HDL_Meta_analysis)


# Function to concatenate duplicate rows based on specified columns
concatenate_duplicates_Chol <- function(data) {
  data %>%
    group_by(Who, Population, StudyID, Sex, HDL_Chol, SerumLevel) %>%
    summarise(across(everything(), ~ {
      # Remove NA values, concatenate unique non-NA values with a comma separator
      unique_vals <- unique(na.omit(.))
      if (length(unique_vals) > 0) {
        paste(unique_vals, collapse = ", ")
      } else {
        NA
      }
    }), .groups = 'drop')
}
# Apply the function to concatenate duplicates
HDL_Meta_analysis <- concatenate_duplicates_Chol(HDL_Meta_analysis)
# Clean and select columns
HDL_Meta_analysis <- dplyr::select(HDL_Meta_analysis, -c(35:37))
dim(HDL_Meta_analysis)
#sum(!is.na(HDL_Meta_analysis$HDL_Chol))
###############################
####################
#########
Population_info_meta_analysis <- import("/Users/Bm0211/RegEx/Water_Energy/Population_Info_Meta_analysis.csv")
colnames(Population_info_meta_analysis)
names(Population_info_meta_analysis)[4]<- "Article_link"

colnames(Population_info_meta_analysis)
dim(Population_info_meta_analysis)
concatenate_duplicates_pop_info <- function(data) {
  data %>%
    group_by(Population, StudyID, Who, Article_link) %>%
    summarise(across(everything(), ~ {
      # Remove NA values, concatenate unique non-NA values with a comma separator
      unique_vals <- unique(na.omit(.))
      if (length(unique_vals) > 0) {
        paste(unique_vals, collapse = ", ")
      } else {
        NA
      }
    }), .groups = 'drop')
}

# Apply the function to concatenate duplicates
Population_info_meta_analysis <- concatenate_duplicates_pop_info(Population_info_meta_analysis)
colnames(Population_info_meta_analysis)
# Clean and select columns
Population_info_meta_analysis <- dplyr::select(Population_info_meta_analysis, c(1:4, 6:11, 17:21, 27, 28, 34))

dim(Population_info_meta_analysis)
###############################
####################
#########
# Standardize AcculturationLevel and Sedentary columns
Population_info_meta_analysis <- Population_info_meta_analysis %>%
  mutate(
    # Standardize AcculturationLevel
    AcculturationLevel = ifelse(str_trim(AcculturationLevel) == "", NA, AcculturationLevel),
    AcculturationLevel = case_when(
      str_detect(str_to_lower(AcculturationLevel), "unacculturated|unacculaturated") ~ "Unacculturated",
      str_detect(str_to_lower(AcculturationLevel), "semi-acculturated|semi_acculturated|semiintegrated|semi-integrated") ~ "Semi-acculturated",
      str_detect(str_to_lower(AcculturationLevel), "isolated") ~ "Isolated",
      str_detect(str_to_lower(AcculturationLevel), "acculturated|integrated|intrgrated") ~ "Integrated",
      str_detect(str_to_lower(AcculturationLevel), "both acculturated and integrated") ~ "Acculturated & Integrated",
      str_detect(str_to_lower(AcculturationLevel), "acculturated") ~ "Acculturated",
      TRUE ~ NA_character_
    ),
    # Standardize Sedentary
    Sedentary = ifelse(str_trim(Sedentary) == "", NA, Sedentary),
    Sedentary = case_when(
      str_detect(str_to_lower(Sedentary), "sedentary/Sedentary") ~ "Sedentary",
      str_detect(str_to_lower(Sedentary), "nomadic|mobile") ~ "Nomadic",
      str_detect(str_to_lower(Sedentary), "semi-nomadic|seminomadic/Semi-nomadic") ~ "Semi-nomadic",
      str_detect(str_to_lower(Sedentary), "active|High activity|Robust physical activity") ~ "Active",
      str_detect(str_to_lower(Sedentary), "Seasonal migration|migration|migrating") ~ "Seasonal Migrant",
      str_detect(str_to_lower(Sedentary), "yes|y") ~ "Yes",
      str_detect(str_to_lower(Sedentary), "no|n") ~ "No",
      TRUE ~ NA_character_
    )
  )
###############################
####################
#########
# Merge datasets
# Function to merge HDL_Meta_analysis with another dataset
merge_with_hdl <- function(hdl_data, other_data, suffix) {
  # Define a helper function to handle column value merging with a separator "|"
  merge_columns <- function(col1, col2) {
    if (is.na(col1) && is.na(col2)) return(NA)
    if (is.na(col1)) return(col2)
    if (is.na(col2)) return(col1)
    if (col1 == col2) return(col1)
    return(paste(col1, col2, sep = "|"))
  }
  # Subset other_data to only include rows with matching StudyID in hdl_data
  other_data <- other_data[other_data$StudyID %in% hdl_data$StudyID, ]
  # Rename columns in other_data to avoid conflict, adding the suffix
  names(other_data) <- ifelse(names(other_data) %in% names(hdl_data) & names(other_data) != "StudyID",
                              paste0(names(other_data), suffix), names(other_data))
  # Perform a left join, keeping only records in HDL_Meta_analysis
  merged_data <- merge(hdl_data, other_data, by = "StudyID", all.x = TRUE)
  # Resolve conflicts for columns that exist in both datasets with the same base name
  common_cols <- intersect(names(hdl_data), gsub(suffix, "", names(other_data)))
  for (col in common_cols) {
    if (col != "StudyID") {
      col_hdl <- paste0(col, "_HDL")
      col_other <- paste0(col, suffix)
      # Check if both columns exist before applying mapply
      if (col_hdl %in% names(merged_data) && col_other %in% names(merged_data)) {
        merged_data[[col]] <- mapply(merge_columns, merged_data[[col_hdl]], merged_data[[col_other]])
        merged_data[[col_hdl]] <- NULL  # Remove the HDL duplicated column
        merged_data[[col_other]] <- NULL  # Remove the other dataset's duplicated column
      }
    }
  }
  
  return(merged_data)
}
# Merge HDL_Meta_analysis with Blood_pressure_meta_analysis
merged_hdl_bp <- merge_with_hdl(HDL_Meta_analysis, Population_info_meta_analysis, "_pop")

dim(merged_hdl_bp)
colnames(merged_hdl_bp)
###############################
####################
#########
# Conversion factors
cholesterol_conversion <- 38.67
triglyceride_conversion <- 88.57

# Columns organized by units columns
cholesterol_groups <- list(
  "SL_Units" = c("SerumLevel", "SL_SD", "SL_SE"),
  "HDL_Units" = c("HDL_Chol", "HDL_SD"),
  "LDL_Units" = c("LDL_Chol", "LDL_SD"),
  "Triglyceride_Units" = c("Triglyceride_Chol", "Triglyceride_SD")
)
# Columns organized by units columns
merged_hdl_bp$SerumLevel <- as.numeric(merged_hdl_bp$SerumLevel)
merged_hdl_bp$SL_SD <- as.numeric(merged_hdl_bp$SL_SD)
merged_hdl_bp$SL_SE <- as.numeric(merged_hdl_bp$SL_SE)
merged_hdl_bp$HDL_Chol <- as.numeric(merged_hdl_bp$HDL_Chol)
merged_hdl_bp$HDL_SD <- as.numeric(merged_hdl_bp$HDL_SD)
merged_hdl_bp$LDL_Chol <- as.numeric(merged_hdl_bp$LDL_Chol)
merged_hdl_bp$LDL_SD <- as.numeric(merged_hdl_bp$LDL_SD)
merged_hdl_bp$Triglyceride_Chol <- as.numeric(merged_hdl_bp$Triglyceride_Chol)
merged_hdl_bp$Triglyceride_SD <- as.numeric(merged_hdl_bp$Triglyceride_SD)

# Conversion function: Convert from mmol/L to mg/dL if unit contains "mmol"
convert_to_mg_dL <- function(value, unit, conversion_factor) {
  if (is.na(value) || is.na(unit)) return(NA)
  unit <- tolower(trimws(unit))  # Standardize by trimming spaces and converting to lowercase
  if (grepl("mmol", unit)) return(value * conversion_factor)  # Matches "mmol" with any spacing
  return(value)  # If not in mmol, return original value
}

# Remove rows with NA in HDL_Chol column
#cleaned_data <- merged_hdl_bp[!is.na(merged_hdl_bp$HDL_Chol), ]
cleaned_data <- merged_hdl_bp
# Iterate over each units column and corresponding data columns
for (units_col in names(cholesterol_groups)) {
  columns <- cholesterol_groups[[units_col]]
  
  # Check if the units column exists in the data
  if (units_col %in% names(cleaned_data)) {
    cleaned_data[[units_col]] <- as.character(cleaned_data[[units_col]])  # Ensure units are character strings
    
    # Loop through each column needing conversion within the group
    for (col in columns) {
      if (col %in% names(cleaned_data)) {
        cleaned_data[[col]] <- as.numeric(as.character(cleaned_data[[col]]))  # Ensure numeric format
        
        # Determine conversion factor based on units column type
        conversion_factor <- if (units_col == "Triglyceride_Units") triglyceride_conversion else cholesterol_conversion
        
        # Create new column name for converted values
        new_col <- paste0(col, "_mg_dL")
        
        # Apply conversion function to each row
        cleaned_data[[new_col]] <- mapply(function(value, unit) 
          convert_to_mg_dL(value, unit, conversion_factor), 
          cleaned_data[[col]], cleaned_data[[units_col]])
      }
    }
  }
}

# Check results
dim(cleaned_data)
# Step 6: Write the cleaned data with old and new columns to a new file
#write.csv(cleaned_data, file = "AYANA_Cleaned_DATA.csv", row.names = FALSE)
# Check final dataset dimensions and structure
colnames(cleaned_data)
dim(cleaned_data)
###############################
####################
#########
##My Meta data - add to the cleaned one above
More_meta <- import("/Users/Bm0211/RegEx/Water_Energy/Lipoproteins_Across_Populations.csv")
colnames(More_meta)
More_meta <- dplyr::select(More_meta, -c(14:27))
names(More_meta)[7] <- "N"
More_meta$N <- as.character(More_meta$N)
More_meta$Lattitude <- as.character(More_meta$Lattitude)
cleaned_data$Lattitude <- as.character(cleaned_data$Lattitude)
str(More_meta)
colnames(More_meta)
dim(More_meta)
# Identify common columns
common_columns <- intersect(colnames(cleaned_data), colnames(More_meta))
# Ensure column names are consistent
colnames(More_meta) <- colnames(More_meta) %>%
  gsub("\n", " ", .) %>% # Replace newline characters with spaces
  trimws() # Trim whitespace

# Combine the datasets with all columns
combined_data <- bind_rows(cleaned_data, More_meta)

# Ensure consistency in column ordering and fill missing values with NA
combined_data <- combined_data %>%
  mutate(across(everything(), ~ ifelse(is.na(.), NA, .)))  # Standardize missing values
dim(combined_data)

##Assigning Ancestry to Populations
##Load file with standardized Ancestry names
Ancestry_names <- import("/Users/Bm0211/RegEx/Water_Energy/AYANA_For_clustering.csv")
head(Ancestry_names)
names(Ancestry_names)[2] <- "Article_access"
# Function to concatenate duplicate rows based on specified columns
concatenate_duplicates_Ancestry <- function(data) {
  data %>%
    group_by(Population, StudyID, Article_access, ANCESTRY) %>%
    summarise(across(everything(), ~ {
      # Remove NA values, concatenate unique non-NA values with a comma separator
      unique_vals <- unique(na.omit(.))
      if (length(unique_vals) > 0) {
        paste(unique_vals, collapse = ", ")
      } else {
        NA
      }
    }), .groups = 'drop')
}
# Apply the function to concatenate duplicates
Ancestry_names <- concatenate_duplicates_Ancestry(Ancestry_names)
names(Ancestry_names)[2]<-"study_ID_Ancestry"
colnames(Ancestry_names)

dim(Ancestry_names)
dim(combined_data)
colnames(combined_data)
# Function to add ancestry information to combined_data based only on Population
add_ancestry <- function(combined_data, ancestry_data) {
  # Remove leading and trailing whitespace from Population column in both datasets
  combined_data <- combined_data %>%
    mutate(Population = trimws(Population))
  ############
  ancestry_data <- ancestry_data %>%
    mutate(Population = trimws(Population)) %>%
    distinct(Population, .keep_all = TRUE)  # Remove duplicates based on Population
  # Perform a left join to add ANCESTRY column based on Population only
  merged_data <- combined_data %>%
    left_join(ancestry_data %>% dplyr::select(Population, ANCESTRY, study_ID_Ancestry), by = "Population")
  # Check if row count matches combined_data; if not, warn about possible mismatches
  if (nrow(merged_data) != nrow(combined_data)) {
    warning("Row count mismatch after merging. Check for unexpected duplicates in ancestry data.")
  }
  return(merged_data)
}
# Add ancestry information
cleaned_data_with_ancestry <- add_ancestry(combined_data, Ancestry_names)

# View the first few rows to verify
dim(cleaned_data_with_ancestry)
colnames(cleaned_data_with_ancestry)
# Define your columns of interest
cols_of_interest <- c(
  "Region", "Lattitude", "N", "SerumLevel", "SL_SD", "Ancestry",
  "HDL_Chol", "HDL_SD", "LDL_Chol", "LDL_SD", "Biome", "ANCESTRY"
)

# 1. Summarize Missingness vs Present
missing_summary <- cleaned_data_with_ancestry %>%
  dplyr::select(all_of(cols_of_interest)) %>%
  summarise(across(everything(), ~sum(is.na(.)))) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Missing") %>%
  mutate(
    Total = nrow(cleaned_data_with_ancestry),
    Present = Total - Missing
  ) %>%
  pivot_longer(cols = c("Missing", "Present"), names_to = "Status", values_to = "Count")

# 2. Plot
missingness_plot <- ggplot(missing_summary, aes(x = reorder(Variable, -Count), y = Count, fill = Status)) +
  geom_col(position = "stack", width = 0.7) +
  geom_text(aes(label = Count), position = position_stack(vjust = 0.5), size = 3, color = "white") +
  scale_fill_manual(values = c("Present" = "forestgreen", "Missing" = "red")) +
  labs(
    title = "Missingness vs Present Data per Variable",
    x = "Variable",
    y = "Number of Entries",
    fill = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 3. Print
print(missingness_plot)

# Save the updated dataset to a CSV file
write.csv(cleaned_data_with_ancestry, file = "/Users/bm0211/RegEx/Water_Energy/2025_Metabolomic_Analysis/Cleaned_data_with_Ancestry_Meta_analysis.csv", row.names = FALSE)

table(cleaned_data_with_ancestry$Population)
cleaned_data_with_ancestry$Population <- make.unique(as.character(cleaned_data_with_ancestry$Population))

#setwd("/Users/bm0211/RegEx/Turkana_CVD")
############################################################################
## Function to Run ANOVA on several groups ##################
FNONCT <- function(x,df1,df2,prob, interval=c(0,10000), my.tol=0.000001){
  temp <- function(ncp) pf(x,df1,df2,ncp) - prob
  return(uniroot(temp, interval, tol = my.tol)$root)
}
############################################################################
power.f <- function(u, n, delta, sig.level=.05){
  fc <- qf(p=sig.level,df1=u,df2=(u+1)*(n-1),lower.tail=FALSE)
  lamda <- (delta^2)*(n*(u+1))
  v <- (u+1)*(n-1)
  
  z1b <- (sqrt(2*(u+lamda)-((u+2*lamda)/(u+lamda)))-
            sqrt((2*v-1)*((u*fc)/v)))/
    sqrt(((u*fc)/v)+((u+2*lamda)/(u+lamda)))
  output <- pnorm(z1b)
  return(output)
}
############################################################################
ind.oneway.second <- function(m, sd, n, unbiased=TRUE, contr=NULL, sig.level=.05, digits=3){
  ##if orthogonal
  if(length(n)==1){
    n <- rep(n, length(m))
  }
  #if biased standard deviation
  if(unbiased==FALSE){
    sd <- ssd2sd(n,sd)
  }
  ############################################################################
  ##(a) anova table
  k  <- length(m)           #number of groups
  Xg  <- sum(n*m)/sum(n)
  dfb <- k - 1              #degree of freedom
  dfw   <- sum(n) - k       #degree of freedom
  MSb <- sum(n * (m - Xg)^2)/(k-1)  #MS between
  MSw <-  sum((n-1)*sd^2)/dfw       #MS within
  SSb <- dfb * MSb
  SSw <- dfw * MSw
  SSt <- SSb + SSw
  f.value <- MSb/MSw                #f value
  anova.table  <- data.frame(matrix(NA,ncol=4, nrow=3))
  rownames(anova.table) <- c("Between (A)", "Within", "Total")
  colnames(anova.table) <- c("SS", "df", "MS", "F")
  anova.table$SS <- c(SSb, SSw,SSt)
  anova.table$df <- c(dfb, dfw, dfb+dfw)
  anova.table$MS <- c(MSb,MSw,NA)
  anova.table$F  <- c(f.value, NA,NA)
  class(anova.table) <- c("anova", "data.frame")
  anova.table <- round(anova.table, digits)
  ############################################################################
  ##(b) omnibus effect size eta and omega squared
  #eta square
  etasq <- SSb / SSt
  delta.lower <- delta.upper <- numeric(length(etasq))
  delta.lower <- try(FNONCT(f.value, dfb, dfw, prob=1-sig.level/2), silent=TRUE)
  delta.upper <- try(FNONCT(f.value, dfb, dfw, prob=sig.level/2), silent=TRUE)
  if(is.character(delta.lower)){
    delta.lower <- 0
  }
  etasq.lower <- delta.lower / (delta.lower + dfb + dfw + 1)
  etasq.upper <- delta.upper / (delta.upper + dfb + dfw + 1)
  ############################################################################
  #omega square
  omegasq <- (SSb - dfb * MSw)/(SSt + MSw)
  sosb_L  <- SSt * etasq.lower
  msw_L   <- (SSt - sosb_L)/dfw
  omegasq.lower <- (sosb_L - (dfb*msw_L))/(SSt+msw_L)
  sosb_U  <- SSt * etasq.upper
  msw_U   <- (SSt - sosb_U)/dfw
  omegasq.upper <- (sosb_U - (dfb*msw_U))/(SSt+msw_U)
  omnibus.es <- round(c(etasq=etasq, etasq.lower=etasq.lower, etasq.upper=etasq.upper),
                      digits)
  ############################################################################
  ##(c) raw contrasts
  temp  <- combinations(k,2)
  cont1 <- matrix(0, nrow=nrow(temp),ncol=k)
  cont1.lab <- rep(0,nrow(temp))
  #in case did not specify contrasts
  for(i in 1:nrow(temp)){
    cont1[i, temp[i,1]] <- 1
    cont1[i, temp[i,2]] <- -1
    cont1.lab[i] <- paste(temp[i,1],"-",temp[i,2], sep="")
    rownames(cont1) <- cont1.lab
  }
  #in case specify contrasts
  if(!is.null(contr)){
    if(is.vector(contr)){
      cont1 <- t(as.matrix(contr))
    }else{
      cont1 <- contr
    }
  }
  ############################################################################
  #F test for contrasts
  psi <-  colSums(t(cont1)  * as.vector(m))                #raw contrasts
  SSpsi <- (psi^2)/colSums(t(cont1^2) / as.vector(n))
  nmat <- matrix(n, nrow=nrow(cont1), ncol=length(n), byrow=TRUE)
  psi.std <- sqrt(MSw * rowSums(cont1 ^ 2/nmat))
  psi.lower <- psi + psi.std * qt(sig.level/2, dfw)
  psi.upper <- psi + psi.std * qt(sig.level/2, dfw, lower.tail=FALSE)
  raw.contrasts <- round(data.frame(mean.diff=psi, lower=psi.lower, upper=psi.upper, std=psi.std), digits)
  rownames(raw.contrasts) <- rownames(cont1)
  ############################################################################
  ##(d) standardized contrasts
  gpsi <- psi/sqrt(MSw)       #effect size
  gpsi.std <- sqrt(rowSums(cont1 ^ 2/nmat))
  gpsi.lower <- gpsi + gpsi.std * qt(sig.level/2, dfw)
  gpsi.upper <- gpsi + gpsi.std * qt(sig.level/2, dfw, lower.tail=FALSE)
  standardized.contrasts <- round(data.frame(es=gpsi, lower=gpsi.lower, upper=gpsi.upper, std=gpsi.std), digits)
  rownames(standardized.contrasts) <- rownames(cont1)
  ##(e) statistical power
  c.delta <- c(.10, .25, .4)
  criterion.power <- round(power.f(sig.level=sig.level, u=dfb, n=sum(n)/k,delta=c.delta), digits)
  names(criterion.power) <- c("small", "medium", "large")
  ##(e) output
  output <- list(anova.table=anova.table, omnibus.es=omnibus.es, raw.contrasts=raw.contrasts, standardized.contrasts = standardized.contrasts, power=criterion.power)
  return(output)
}
###############################################################################################
##Fake Data to test the Function
mean <- c(90,85,92,100,102,106)
sd <- c(9.035613,11.479667,9.760268,7.662572,9.830258,9.111457)
no_of_participants <- c(9,9,9,9,9,9)
ind.oneway.second (mean, sd, no_of_participants) ##This this the function and compares Max of 7 groups at once
###############################################################################################
# Load CSV file and convert columns to Numeric then loop through all the populations
df <- import("/Users/bm0211/RegEx/Water_Energy/2025_Metabolomic_Analysis/Meta_analysis_Frozen_April_2025.csv")
# Define biome standardization mapping
df <- df %>%
  mutate(Biome = str_trim(Biome)) %>%
  mutate(Biome = ifelse(Biome == "" | Biome == " ", NA, Biome))
biome_cleaning_dict <- c(
  "temperate rainforest" = "forest",
  "Temperate Rainforest" = "forest",
  "Temperate rainforest" = "forest",
  "temperate rainfoest" = "forest",
  "tropical rainforest" = "forest",
  `tropical rainforest?` = "forest",
  "Tropical rainforest" = "forest",
  "Tropical Rainforest" = "forest",
  "tropical and subtropical moist broadleaf forests" = "forest",
  "subtripical rainforest" = "forest",
  "subtropical forest" = "forest",
  "temperate deciduous forest" = "Temperate Forest",
  "Temperate forest" = "forest",
  "temperate forest" = "forest",
  "Swamp forest" = "Wetland",
  "taiga" = "Taiga",
  "Temperate Shrubland" = "Temperate Shrubland",
  "temperate shrubland" = "Temperate Shrubland",
  `temperate shrubland/forest` = "Temperate Shrubland",
  `Temperate shrubland/forest` = "Temperate Shrubland",
  "Grassland" = "Temperate Grassland",
  "grassland" = "Temperate Grassland",
  "grasslands" = "Temperate Grassland",
  "Temperate Grassland" = "Temperate Grassland",
  "Desert" = "Desert",
  "tropical desert" = "Desert",
  `Warm-temperate desert` = "Desert",
  "urban" = "Urban",
  "Arctic Tundra" = "Tundra",
  "arctic tundra" = "Tundra",
  "Tundra" = "Tundra",
  "tundra" = "Tundra",
  `alpine/tundra` = "Tundra",
  `Arctic-Alpine desert` = "Tundra",
  "Savanna" = "Savanna",
  "savanna" = "Savanna",
  "savannah" = "Savanna",
  "Woodland" = "Woodland",
  "woodland" = "Woodland",
  `, woodland` = "Woodland",
  "tropical" = "Tropical", 
  "temperate forest" = "forest",
  "tropical rainforest" = "forest",# too vague, drop or reassign manually
  "NA" = NA,
  " " = NA
)
# Apply cleaning to your dataframe
df <- df %>%
  mutate(
    Biome = recode(Biome, !!!biome_cleaning_dict),
    Biome = str_trim(tolower(Biome)), 
    Biome = factor(Biome)  # Optional: convert to factor
  )
table(df$Biome)
# Clean Region column and standardize names
df <- df %>%
  mutate(
    Region = case_when(
      Region %in% c("north america", "north american", "north america") ~ "North America",
      Region %in% c("south america", "sourth america") ~ "South America",
      Region %in% c("asia pacific region") ~ "Asia Pacific",
      Region %in% c("west africa") ~ "West Africa",
      Region == "central america" ~ "Central America",
      Region == "middle east" ~ "Middle East",
      Region == "distinct" ~ "Distinct",
      Region == "africa" ~ "Africa",
      Region == "asia" ~ "Asia",
      Region == "europe" ~ "Europe",
      Region == "australia" ~ "Australia",
      Region == "pacific" ~ "Pacific",
      Region == "oceania" ~ "Oceania",
      Region == "na" | Region == "" ~ NA_character_,
      TRUE ~ stringr::str_to_title(Region)  # fallback, capitalizes words
    ),
    Region = tolower(trimws(Region)),
  )
##
##
df <- df %>%
  mutate(
    Latitude_Band = case_when(
      is.na(Lattitude) ~ "Unknown",
      Lattitude < -15 ~ "Far South",
      Lattitude >= -15 & Lattitude < -10 ~ "South",
      Lattitude >= -10 & Lattitude <= 10 ~ "Equatorial",
      Lattitude > 10 & Lattitude <= 30 ~ "North",
      Lattitude > 30 ~ "Far North"
    )
  )
#####################################################################
colnames(df)
df$Lattitude <- as.numeric(df$Lattitude)
df$N <- as.numeric(df$N)
df$LDL_SD_mg_dL <- as.numeric(df$LDL_SD_mg_dL)
df$LDL_Chol_mg_dL <- as.numeric(df$LDL_Chol_mg_dL)
df$HDL_Chol_mg_dL <- as.numeric(df$HDL_Chol_mg_dL)
df$HDL_SD_mg_dL <- as.numeric(df$HDL_SD_mg_dL)
df$SerumLevel_mg_dL <- as.numeric(df$SerumLevel_mg_dL)
df$SL_SD_mg_dL <- as.numeric(df$SL_SD_mg_dL)
df$Population <- make.unique(as.character(df$Population))
## Compare all to Framinghal
df$Population[df$Population == "Framingham study"] <- "Framingham_study"

ref_name <- "Framingham_study" #or "TURKANA_pastoralist"or Framingham study
str(df)
table(df$Population)
# Function to run batches
run_anova_batched <- function(df, biomarker_mean, biomarker_sd, n_col, batch_size = 9) {
  df_clean <- df %>%
    dplyr::select(Population, N = {{n_col}}, Mean = {{biomarker_mean}}, SD = {{biomarker_sd}}) %>%
    filter(!is.na(Mean) & !is.na(SD) & !is.na(N)) %>%
    #filter(Mean > 10 & Mean <= 230 & SD > 7) %>%
    arrange(desc(Population == ref_name))  # Ensure reference is first
  ref_row <- df_clean %>% filter(Population == ref_name)
  other_rows <- df_clean %>% filter(Population != ref_name)
  batches <- split(other_rows, (seq_len(nrow(other_rows)) - 1) %/% batch_size)
  # Safe version of the function
  safe_ind_anova <- purrr::safely(function(batch, i) {
    dat <- bind_rows(ref_row, batch)
    m <- dat$Mean
    sd <- dat$SD
    n <- dat$N
    popnames <<- dat$Population  # needed for extract_vs_reference()
    cat("\n🧪 Batch index:", i, "\n")
    print(data.frame(Population = popnames, Mean = m, SD = sd, N = n))
    ind.oneway.second(m, sd, n)
  })
  # Apply the safe function
  results_raw <- map2(batches, seq_along(batches), safe_ind_anova)
  # Filter to only keep successful results
  results <- map(results_raw, "result") %>% compact()
  return(results)
}
# Function to extract standardized comparisons
extract_vs_reference <- function(result_batch) {
  std_contrasts <- result_batch$standardized.contrasts
  ref_contrasts <- std_contrasts[grepl("^1-", rownames(std_contrasts)), ]
  idx <- as.numeric(sub("1-", "", rownames(ref_contrasts)))
  data.frame(
    Population = popnames[idx],
    SMD = ref_contrasts$es,
    CI_lower = ref_contrasts$lower,
    CI_upper = ref_contrasts$upper
  )
}
# Helper to compute significance
extract_all_smd <- function(results, df, mean_col, sd_col) {
  map2_dfr(results, seq_along(results), function(result, i) {
    batch_pops <- df %>%
      filter(!is.na({{mean_col}}) & !is.na({{sd_col}})) %>%
      arrange(desc(Population == ref_name)) %>%
      slice(((i - 1) * 9 + 1):(i * 9)) %>%
      pull(Population)
    popnames <<- c("Framingham_study", make.unique(as.character(batch_pops)))
    extract_vs_reference(result)
  }) %>%
    mutate(
      significant = ifelse((SMD > 1 & CI_lower > 1) | (SMD < -1 & CI_upper < -1), TRUE, FALSE)
    )
}
# Run comparisons
colnames(df)
hdl_results <- run_anova_batched(df, HDL_Chol_mg_dL, HDL_SD_mg_dL, N)
ldl_results <- run_anova_batched(df, LDL_Chol_mg_dL, LDL_SD_mg_dL, N)
chol_results <- run_anova_batched(df, SerumLevel_mg_dL, SL_SD_mg_dL, N)
############################################################################
hdl_smd <- extract_all_smd(hdl_results, df, HDL_Chol_mg_dL, HDL_SD_mg_dL)
ldl_smd <- extract_all_smd(ldl_results, df, LDL_Chol_mg_dL, LDL_SD_mg_dL)
chol_smd <- extract_all_smd(chol_results, df, SerumLevel_mg_dL, SL_SD_mg_dL)
############################################################################
highlight_pops <- c("UK Biobank - Healthy", "TURKANA_pastoralist")
############################################################################
plot_smd <- function(data, title, sig_color, invert_order = FALSE) {
  data <- data %>% filter(SMD >= -5 & SMD <= 5)
  ordering <- if (invert_order) -data$SMD else data$SMD
  data <- data %>% mutate(Pop_ordered = reorder(Population, ordering))
  ggplot(data, aes(x = SMD, y = Pop_ordered)) +
    geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper, color = significant), height = 0.25) +
    geom_point(aes(color = significant), size = 1) +
    scale_color_manual(values = c("TRUE" = sig_color, "FALSE" = "black"), guide = "none") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(-1, 1), linetype = "solid", color = sig_color) +
    labs(title = title, x = "Cohen's d", y = "Population") +
    ggrepel::geom_text_repel(
      data = data %>% filter(Population %in% highlight_pops),
      aes(label = Population, color = significant),
      size = 3.5,
      nudge_x = 0.3,
      direction = "y",
      hjust = 0,
      segment.color = "grey50",
      max.overlaps = Inf
    ) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank(),
      axis.text.y = element_text(size = 3),
      plot.title = element_text(face = "bold", size = 10)
    )
}
############################################################################
# Pie chart function
plot_pie <- function(data, title, color) {
  sig_summary <- data %>%
    mutate(Significance = ifelse(significant, "Significant", "Not Significant")) %>%
    count(Significance) %>%
    mutate(perc = round(n / sum(n) * 100, 1),
           label = paste0(perc, "%"))
  
  ggplot(sig_summary, aes(x = "", y = n, fill = Significance)) +
    geom_col(width = 1, color = "white") +
    coord_polar(theta = "y") +
    geom_text(aes(label = label), position = position_stack(vjust = 0.5), color = "white", size = 5) +
    scale_fill_manual(values = c("Significant" = color, "Not Significant" = "black")) +
    labs(title = title, fill = NULL) +
    theme_void(base_size = 13) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}
# Generate plots with adjusted ordering
p1 <- plot_smd(hdl_smd, "HDL All vs. Framingham", "#A7BA42", invert_order = TRUE)
p2 <- plot_smd(ldl_smd, "LDL All vs. Framingham", "#DFC41B", invert_order = TRUE)
p3 <- plot_smd(chol_smd, "Chol All vs. Framingham", "red", invert_order = TRUE)
# Combine
combined_smd <- p1 + p2 + p3
print(combined_smd)
ggsave("/Users/bm0211/RegEx/Water_Energy/2025_Metabolomic_Analysis/SMD_vs_TURKANA_pastoralist.png", combined_smd, width = 16, height = 16, dpi = 600)

combined_pie <- plot_pie(hdl_smd, "HDL Significance", "#A7BA42") +
  plot_pie(ldl_smd, "LDL Significance", "#DFC41B") +
  plot_pie(chol_smd, "Chol Significance", "red")
print(combined_pie)
ggsave("/Users/bm0211/RegEx/Water_Energy/2025_Metabolomic_Analysis/Significance_PieCharts.png", combined_pie, width = 12, height = 6, dpi = 600)
###
###
############################################################################
# Clean and filter to valid entries
mean_distributions <- df %>%
  filter(!is.na(LDL_Chol_mg_dL) & !is.na(HDL_Chol_mg_dL) & !is.na(SerumLevel_mg_dL))
# Plot for HDL
p_hdl <- ggplot(mean_distributions, aes(x = HDL_Chol_mg_dL)) +
  geom_histogram(binwidth = 2, fill = "#A7BA42", color = "white", alpha = 0.9) +
  labs(
    title = "Distribution of HDL Means Across Populations",
    x = "HDL (mg/dL)", y = "Number of Populations"
  ) +
  theme_minimal(base_size = 13)
# Plot for LDL
p_ldl <- ggplot(mean_distributions, aes(x = LDL_Chol_mg_dL)) +
  geom_histogram(binwidth = 5, fill = "#DFC41B", color = "white", alpha = 0.9) +
  labs(
    title = "Distribution of LDL Means Across Populations",
    x = "LDL (mg/dL)", y = "Number of Populations"
  ) +
  theme_minimal(base_size = 13)
# Plot for Cholesterol
p_chol <- ggplot(mean_distributions, aes(x = SerumLevel_mg_dL)) +
  geom_histogram(binwidth = 5, fill = "red", color = "white", alpha = 0.9) +
  labs(
    title = "Distribution of Total Cholesterol Means",
    x = "Cholesterol (mg/dL)", y = "Number of Populations"
  ) +
  theme_minimal(base_size = 13)
# Combine plots
combined_dist <- p_hdl / p_ldl / p_chol  # stack vertically
print(combined_dist)
# Save output
ggsave(
  "/Users/bm0211/RegEx/Water_Energy/2025_Metabolomic_Analysis/Distributions_of_Means.png",
  combined_dist,
  width = 10, height = 12, dpi = 600
)
########################################################################################################################
colnames(df)
# Add scaled lipid traits
df <- df %>%
  mutate(
    HDL_scaled = rescale(HDL_Chol_mg_dL, to = c(0, 1), na.rm = TRUE),
    LDL_scaled = rescale(LDL_Chol_mg_dL, to = c(0, 1), na.rm = TRUE),
    Chol_scaled = rescale(SerumLevel_mg_dL, to = c(0, 1), na.rm = TRUE)
  )
# Function to remove outliers based on IQR
remove_outliers <- function(df, var) {
  q1 <- quantile(df[[var]], 0.25, na.rm = TRUE)
  q3 <- quantile(df[[var]], 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  lower <- q1 - 1.5 * iqr
  upper <- q3 + 1.5 * iqr
  df[df[[var]] >= lower & df[[var]] <= upper, ]
}
# Custom color map for each trait
trait_colors <- c(
  "HDL" = "#33a02c",         # Green
  "LDL" = "#ffcf00",         # Yellow
  "Cholesterol" = "#e31a1c"  # Red
)
df <- df %>% mutate(abs_Lattitude = abs(Lattitude))
# Plotting function with curved fit and original (signed) latitude
plot_scaled_vs_latitude <- function(df, scaled_var, label) {
  df_clean <- remove_outliers(df, scaled_var)
  formula <- as.formula(paste(scaled_var, "~ abs_Lattitude"))
  model <- lm(formula, data = df_clean)
  summary_model <- summary(model)
  r2 <- round(summary_model$r.squared, 3)
  pval <- signif(coef(summary_model)[2, 4], 3)
  color <- trait_colors[label]
  ggplot(df_clean, aes_string(x = "abs_Lattitude", y = scaled_var)) +
    geom_point(alpha = 0.6, color = color) +
    geom_smooth(method = "loess", se = TRUE, color = "black", span = 1) +
    labs(
      title = paste(label, "vs Latitude"),
      subtitle = paste("R² =", r2, "| p =", pval),
      x = "Latitude (°) Absolute value",
      y = paste("Scaled", label, "(0–1)")
    ) +
    theme_minimal(base_size = 14)
}
# Generate individual plots
p1 <- plot_scaled_vs_latitude(df, "HDL_scaled", "HDL")
p2 <- plot_scaled_vs_latitude(df, "LDL_scaled", "LDL")
p3 <- plot_scaled_vs_latitude(df, "Chol_scaled", "Cholesterol")
# Combine side-by-side
combined_plot22 <- p1 + p2 + p3 + plot_layout(ncol = 3)
print(combined_plot22)
# Save
ggsave(
  "/Users/bm0211/RegEx/Water_Energy/2025_Metabolomic_Analysis/Scaled_Lipids_vs_Latitude.png",combined_plot22,width = 16, height = 5, dpi = 600)

########################################################################################################################
########################################################################################################################
### Heritability Meta_analysis
### Heritability Meta_analysis
### Heritability Meta_analysis
### Heritability Meta_analysis
########################################################################################################################
# Turkana data from GCTA results.............. GCTA_Heritability.sh 
turkana <- data.frame(
  Population = "Turkana",
  Trait = c("Chol", "LDL", "HDL"),
  h2 = c(0.477520, 0.185399, 0.413902),
  se = c(0.256083, 0.315208, 0.223360),
  N = c(378, 296, 416))
########################################################################################################################
### Heritability Meta_analysis
# Turkana data from newest EMMA.vcf shared by kristina on slack
turkana <- data.frame(
  Population = "Turkana",
  Trait = c("Chol", "LDL", "HDL"),
  h2 = c(0.75216, 0.843867, 0.232369),
  se = c(0.172902, 0.21615, 0.204967),
  N = c(747, 664, 768))
########################################################################################################################
# Momin et al. 2023 data
########################################################################################################################
momin_data <- data.frame(
  Population = rep(c("White British", "Other European", "South Asian", "African"), each = 3),
  Trait = rep(c("Chol", "LDL", "HDL"), times = 4),
  h2 = c(0.256, 0.241, 0.310,
         0.278, 0.260, 0.308,
         0.162, 0.155, 0.317,
         0.308, 0.291, 0.287),
  se = c(0.012, 0.011, 0.013,
         0.016, 0.015, 0.017,
         0.035, 0.038, 0.042,
         0.035, 0.037, 0.038),
  N = c(30000, 30000, 30000,
        26457, 26457, 26457,
        6199, 6199, 6199,
        6179, 6179, 6179)
)
#####################################
# Hispanic/Latino (HCHS/SOL
# New data: Latin American population from screenshot
# CI = [lower, upper], SE estimated as (upper - lower)/3.92 (approximate for 95% CI)
# https://doi.org/10.1186/s12944-017-0591-6
#####################################
latin_data <- data.frame(
  Population = "Latin America",
  Trait = c("Chol", "LDL", "HDL"),
  h2 = c(0.22, 0.21, 0.24),
  se = c( (0.29 - 0.15) / 3.92,
          (0.27 - 0.14) / 3.92,
          (0.30 - 0.17) / 3.92),
  N = c(10264, 10264, 10264))
#####################################
# Brazil
## 10.1186/1471-2350-9-32 
#####################################
brazil_data <- data.frame(
  Population = "Brazil",
  Trait = c("Chol", "LDL", "HDL"),
  h2 = c(0.29, 0.26, 0.32),
  se = c(0.049, 0.049, 0.046),
  N = c(1666, 1666, 1666))
#####################################
# China Twins
# https://doi.org/10.3390/nu15010164
china_twin_data <- data.frame(
  Population = "Chinese twins",
  Trait = c("Chol", "LDL", "HDL"),
  h2 = c(0.63, 0.60, 0.69),
  se = c( (0.67 - 0.59) / 3.92, 
                  (0.64 - 0.55) / 3.92, 
                  (0.72 - 0.65) / 3.92),
  N = c(2842, 2842, 2842))
####################################
# us_veteran_data
# https://doi.org/10.1375/twin.10.5.703
us_veteran_data <- data.frame(
  Population = "US_Veteran_Male_Twins",
  Trait      = c("Chol", "LDL", "HDL"),
  h2         = c(0.46,   0.49,  0.50),  # Exam 1 heritabilities
  se         = c(NA,     NA,    NA),    # not reported in paper
  N          = c(1025,   1025,  1025))
####################################
# Princeton Ohio Follow Up
# https://www.jlr.org/article/S0022-2275(20)35256-1/fulltext
pfs_data <- data.frame(
  Population = rep(c("USA_White_LRC",
                     "USA_White_PFS",
                     "USA_AA_LRC",
                     "USA_AA_PFS"), each = 3),
  Trait = rep(c("Chol", "LDL", "HDL"), times = 4),
  h2 = c(
    0.70, 0.65, 0.75,  # White, LRC (1970s)
    0.54, 0.46, 0.44,  # White, PFS (2000s)
    0.22, 0.17, 0.61,  # African American, LRC
    0.25, 0.27, 0.48   # African American, PFS
  ),
  se = c(
    0.07, 0.07, 0.07,  # White, LRC
    0.08, 0.08, 0.07,  # White, PFS
    0.15, 0.17, 0.16,  # African American, LRC
    0.14, 0.13, 0.13   # African American, PFS
  ),
  N = c(
    rep(1077, 6),  # white sample size (parents + offspring)
    rep(367, 6)    # African American sample size
  )
)

####################################
# Busselton Family Heart Study - European Ancestry
# https://www.jlr.org/article/S0022-2275(20)43503-5/fulltext
busselton_data <- data.frame(
  Population = "Busselton Study - European Ancestry",
  Trait      = c("Chol", "LDL", "HDL"),
  h2         = c(0.57, 0.52, 0.59),
  se         = c(0.03, 0.04, 0.03),
  N          = c(4492, 4492, 4492)
)

# Combine all data
all_data <- bind_rows(turkana, momin_data, latin_data, brazil_data, china_twin_data, us_veteran_data, pfs_data, busselton_data)

# Order traits
all_data$Trait <- factor(all_data$Trait, levels = c("HDL", "LDL", "Chol", "Height"))

# Custom population colors
#pop_colors <- c(
 # "White British" = "#C98B9A",
  #"Turkana" = "#8B0000",
  #"South Asian" = "cyan",
  #"Other European" = "black",
  #"African" = "gray70",
  #"Latin American" = "#009E73"
#)
####################################
# Plot
# Create color vector: Turkana highlighted, everyone else black
####################################
pop_levels <- unique(all_data$Population)
pop_colors <- setNames(
  ifelse(pop_levels == "Turkana", "red", "grey"),
  pop_levels
)

ggplot(all_data, aes(x = Trait, y = h2, fill = Population)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.9),width = 0.72) +  
  geom_errorbar(aes(ymin = h2 - se, ymax = h2 + se),
                width = 0.2,
                position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = pop_colors) +
  labs(title = "Heritability Across Global Cohorts",
       x = "Trait",
       y = "Heritability (h²)",
       fill = "Cohorts included") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right") + coord_flip()

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
# 1. HDL Cholesterol Plot ----------------------------------------------
hdl_data <- cleaned_data_with_ancestry %>%
  filter(!is.na(HDL_Chol_mg_dL))
colnames(hdl_data)
hdl_plot <- ggplot(hdl_data, 
                   aes(x = HDL_Chol_mg_dL, 
                       y = reorder(Population, HDL_Chol_mg_dL, FUN = base::mean))) +
  # Default error bars for HDL 
  geom_errorbar(aes(xmin = HDL_Chol_mg_dL - HDL_SD_mg_dL, 
                    xmax = HDL_Chol_mg_dL + HDL_SD_mg_dL),
                width = 0, color = "#A7BA42", 
                position = position_jitter(height = 0.1)) +
  # Default mean points for HDL
  geom_point(size = 1, color = "#A7BA42",
             position = position_jitter(height = 0.1)) +
  # Overplot TURKANA error bars in wine red
  geom_errorbar(data = subset(hdl_data, Population == "TURKANA_pastoralist"),
                aes(xmin = HDL_Chol_mg_dL - HDL_SD_mg_dL, 
                    xmax = HDL_Chol_mg_dL + HDL_SD_mg_dL),
                width = 0, color = "#722F37", size = 2,
                position = position_jitter(height = 0.1)) +
  # Overplot TURKANA mean points in wine
  geom_point(data = subset(hdl_data, Population == "TURKANA_pastoralist"),
             size = 5, color = "#722F37",
             position = position_jitter(height = 0.1)) +
  # Overplot framingham error bars in blue
  geom_errorbar(data = subset(hdl_data, Population == "Framingham study"),
                aes(xmin = HDL_Chol_mg_dL - HDL_SD_mg_dL, 
                    xmax = HDL_Chol_mg_dL + HDL_SD_mg_dL),
                width = 0, color = "blue", size = 2,
                position = position_jitter(height = 0.1)) +
  # Overplot framingham mean points in blue
  geom_point(data = subset(hdl_data, Population == "Framingham study"),
             size = 5, color = "blue",
             position = position_jitter(height = 0.1)) +
  # Overplot UK Biobank - Healthy error bars in pink
  geom_errorbar(data = subset(hdl_data, Population == "UK Biobank - Healthy"),
                aes(xmin = HDL_Chol_mg_dL - HDL_SD_mg_dL, 
                    xmax = HDL_Chol_mg_dL + HDL_SD_mg_dL),
                width = 0, color = "#C98B9A", size = 2,
                position = position_jitter(height = 0.1)) +
  # Overplot UK Biobank - Healthy mean points in pink
  geom_point(data = subset(hdl_data, Population == "UK Biobank - Healthy"),
             size = 5, color = "#C98B9A",
             position = position_jitter(height = 0.1)) +
  labs(x = "HDL Cholesterol (mg/dL)", y = "Population",
       title = "HDL Cholesterol by Population") +
  theme_minimal() +
  theme(panel.grid = element_blank())

## Print
print(hdl_plot)
# 2. LDL Cholesterol Plot ----------------------------------------------
ldl_data <- cleaned_data_with_ancestry %>%
  filter(!is.na(LDL_Chol_mg_dL))

ldl_plot <- ggplot(ldl_data, 
                   aes(x = LDL_Chol_mg_dL, 
                       y = reorder(Population, LDL_Chol_mg_dL, FUN = base::mean))) +
  # Default error bars for LDL (yellowish)
  geom_errorbar(aes(xmin = LDL_Chol_mg_dL - LDL_SD_mg_dL, 
                    xmax = LDL_Chol_mg_dL + LDL_SD_mg_dL),
                width = 0, color = "#DFC41B", 
                position = position_jitter(height = 0.1)) +
  # Default mean points for LDL (yellowish)
  geom_point(size = 1, color = "#DFC41B",
             position = position_jitter(height = 0.1)) +
  # Overplot TURKANA error bars in wine
  geom_errorbar(data = subset(ldl_data, Population == "TURKANA_pastoralist"),
                aes(xmin = LDL_Chol_mg_dL - LDL_SD_mg_dL, 
                    xmax = LDL_Chol_mg_dL + LDL_SD_mg_dL),
                width = 0, color = "#722F37", size = 2,
                position = position_jitter(height = 0.1)) +
  # Overplot TURKANA mean points in wine
  geom_point(data = subset(ldl_data, Population == "TURKANA_pastoralist"),
             size = 5, color = "#722F37",
             position = position_jitter(height = 0.1)) +
  # Overplot framingham error bars in blue
  geom_errorbar(data = subset(ldl_data, Population == "Framingham study"),
                aes(xmin = LDL_Chol_mg_dL - LDL_SD_mg_dL, 
                    xmax = LDL_Chol_mg_dL + LDL_SD_mg_dL),
                width = 0, color = "blue", size = 2,
                position = position_jitter(height = 0.1)) +
  # Overplot framingham mean points in blue
  geom_point(data = subset(ldl_data, Population == "Framingham study"),
             size = 5, color = "blue",
             position = position_jitter(height = 0.1)) +
  # Overplot UK Biobank - Healthy error bars in pink
  geom_errorbar(data = subset(ldl_data, Population == "UK Biobank - Healthy"),
                aes(xmin = LDL_Chol_mg_dL - LDL_SD_mg_dL, 
                    xmax = LDL_Chol_mg_dL + LDL_SD_mg_dL),
                width = 0, color = "#C98B9A", size = 2,
                position = position_jitter(height = 0.1)) +
  # Overplot UK Biobank - Healthy mean points in pink
  geom_point(data = subset(ldl_data, Population == "UK Biobank - Healthy"),
             size = 5, color = "#C98B9A",
             position = position_jitter(height = 0.1)) +
  labs(x = "LDL Cholesterol (mg/dL)", y = "Population",
       title = "LDL Cholesterol by Population") +
  theme_minimal() +
  theme(panel.grid = element_blank())

## Print
print(ldl_plot)
# 3. Total Cholesterol Plot ----------------------------------------------

chol_data <- cleaned_data_with_ancestry %>%
  filter(!is.na(SerumLevel_mg_dL))

chol_plot <- ggplot(chol_data, 
                    aes(x = SerumLevel_mg_dL, 
                        y = reorder(Population, SerumLevel_mg_dL, FUN = base::mean))) +
  # Default error bars for total Cholesterol (black)
  geom_errorbar(aes(xmin = SerumLevel_mg_dL - SL_SD_mg_dL, 
                    xmax = SerumLevel_mg_dL + SL_SD_mg_dL),
                width = 0, color = "black", 
                position = position_jitter(height = 0.1)) +
  # Default mean points for total Cholesterol (black)
  geom_point(size = 1, color = "black",
             position = position_jitter(height = 0.1)) +
  # Overplot TURKANA error bars in wine
  geom_errorbar(data = subset(chol_data, Population == "TURKANA_pastoralist"),
                aes(xmin = SerumLevel_mg_dL - SL_SD_mg_dL, 
                    xmax = SerumLevel_mg_dL + SL_SD_mg_dL),
                width = 0, color = "#722F37", size = 2,
                position = position_jitter(height = 0.1)) +
  # Overplot TURKANA mean points in wine
  geom_point(data = subset(chol_data, Population == "TURKANA_pastoralist"),
             size = 5, color = "#722F37",
             position = position_jitter(height = 0.1)) +
  # Overplot framingham error bars in blue
  geom_errorbar(data = subset(chol_data, Population == "Framingham study"),
                aes(xmin = SerumLevel_mg_dL - SL_SD_mg_dL, 
                    xmax = SerumLevel_mg_dL + SL_SD_mg_dL),
                width = 0, color = "blue", size = 2,
                position = position_jitter(height = 0.1)) +
  # Overplot framingham mean points in blue
  geom_point(data = subset(chol_data, Population == "Framingham study"),
             size = 5, color = "blue",
             position = position_jitter(height = 0.1)) +
  # Overplot UK Biobank - Healthy error bars in pink
  geom_errorbar(data = subset(chol_data, Population == "UK Biobank - Healthy"),
                aes(xmin = SerumLevel_mg_dL - SL_SD_mg_dL, 
                    xmax = SerumLevel_mg_dL + SL_SD_mg_dL),
                width = 0, color = "#C98B9A", size = 2,
                position = position_jitter(height = 0.1)) +
  # Overplot UK Biobank - Healthy mean points in pink
  geom_point(data = subset(chol_data, Population == "UK Biobank - Healthy"),
             size = 5, color = "#C98B9A",
             position = position_jitter(height = 0.1)) +
  labs(x = "Total Cholesterol (mg/dL)", y = "Population",
       title = "Total Cholesterol by Population") +
  theme_minimal() +
  theme(panel.grid = element_blank())
print(chol_plot)
# Save the plots as PNG files:

ggsave("HDL_Cholesterol_by_Population.png", plot = hdl_plot, width = 10, height = 8, dpi = 300)
ggsave("LDL_Cholesterol_by_Population.png", plot = ldl_plot, width = 10, height = 8, dpi = 300)
ggsave("Total_Cholesterol_by_Population.png", plot = chol_plot, width = 10, height = 8, dpi = 300)


# Read the cleaned data with ancestry information after manual check
Cleaned_datamanually_checked <- import("Cleaned_data_with_Ancestry_manually_checked.csv")
## Cleaned_datamanually_checked <- cleaned_data_with_ancestry #The Data I wrote for Manual Cleaning
dim(Cleaned_datamanually_checked)
table(Cleaned_datamanually_checked$Population)
colnames(Cleaned_datamanually_checked)
# Clean and process the data
# Convert the four columns to numeric, coercing non-numeric values to NA
Cleaned_datamanually_checked <- Cleaned_datamanually_checked %>%
  mutate(
    N = as.numeric(N),
    SerumLevel_mg_dL = as.numeric(SerumLevel_mg_dL),
    HDL_Chol_mg_dL = as.numeric(HDL_Chol_mg_dL),
    LDL_Chol_mg_dL = as.numeric(LDL_Chol_mg_dL)
  )
# Filter out rows where any of the three columns of interest is NA
Cleaned_datamanually_checked <- Cleaned_datamanually_checked %>%
  filter(!is.na(N) & !is.na(SerumLevel_mg_dL) & !is.na(HDL_Chol_mg_dL) & !is.na(LDL_Chol_mg_dL))
# Function to calculate pooled mean
pooled_mean <- function(values, weights) {
  sum(values * weights, na.rm = TRUE) / sum(weights, na.rm = TRUE)
}
# Clean and process the data
Cleaned_datamanually_checked <- Cleaned_datamanually_checked %>%
  group_by(study_ID_Ancestry) %>%
  summarize(
    # Calculate pooled means for each biomarker
    SerumLevel_mg_dL = pooled_mean(SerumLevel_mg_dL, N),
    HDL_Chol_mg_dL = pooled_mean(HDL_Chol_mg_dL, N),
    LDL_Chol_mg_dL = pooled_mean(LDL_Chol_mg_dL, N),
    N = sum(N, na.rm = TRUE),
    # Merge non-biomarker columns
    across(
      everything(),
      ~ {
        if (is.numeric(.)) . <- as.character(.)
        if (n_distinct(na.omit(.)) == 1) unique(na.omit(.)) else paste(unique(na.omit(.)), collapse = " | ")
      },
      .names = "merged_{col}"
    )
  ) %>%
  ungroup() %>%
  dplyr::select(-starts_with("merged_"), everything())
dim(Cleaned_datamanually_checked)
# Save or view the cleaned data
write.csv(Cleaned_datamanually_checked, "Cleaned_datamanually_checked.csv", row.names = FALSE)
print(Cleaned_datamanually_checked)

# Identify rows with complete data for SerumLevel, HDL, and LDL cholesterol
rows_with_complete_data <- Cleaned_datamanually_checked %>%
  dplyr::select(SerumLevel_mg_dL, HDL_Chol_mg_dL, LDL_Chol_mg_dL) %>%
  complete.cases()
# Subset data with complete rows only and select relevant columns
data_for_clustering <- Cleaned_datamanually_checked[rows_with_complete_data, ] %>%
  dplyr::select(SerumLevel_mg_dL, HDL_Chol_mg_dL, LDL_Chol_mg_dL)
# Create a new variable for the LDL/HDL ratio (LDL is the numerator, HDL is the denominator)
data_for_clustering$LDL_HDL_Ratio <- data_for_clustering$LDL_Chol_mg_dL / data_for_clustering$HDL_Chol_mg_dL
# For clustering, use SerumLevel and LDL_HDL_Ratio only
data_for_clustering_reduced <- data_for_clustering %>%
  dplyr::select(SerumLevel_mg_dL, LDL_HDL_Ratio)
# Scale the unweighted data and add the scaled dimensions back into the original dataframe
scaled_unweighted_data <- scale(data_for_clustering_reduced)
clustered_data_unweighted <- Cleaned_datamanually_checked[rows_with_complete_data, ]
clustered_data_unweighted <- cbind(clustered_data_unweighted, 
                                   Dimension1 = scaled_unweighted_data[, 1], 
                                   Dimension2 = scaled_unweighted_data[, 2])
# Run K-means clustering on the scaled unweighted data
set.seed(123)  # For reproducibility
kmeans_result_unweighted <- kmeans(scaled_unweighted_data, centers = 2, nstart = 50)
clustered_data_unweighted$Cluster <- as.factor(kmeans_result_unweighted$cluster)

# Rename the coloring legend column (assuming column 41 holds the Study_Population information)
##Use the Ancestry column as the study population, I used it because there is uniform naming for easy graphing
names(clustered_data_unweighted)[41] <- "Study_Population"
colnames(clustered_data_unweighted)
table(clustered_data_unweighted$merged_ANCESTRY)
# Create the plot object with updated axis labels reflecting the LDL/HDL ratio
plot_unweighted <- ggplot(clustered_data_unweighted, 
                          aes(x = Dimension1, y = Dimension2, color = Study_Population)) +
  geom_point(size = 2) +
  stat_ellipse(aes(group = Cluster), type = "t", level = 0.95, alpha = 0.2, color = "grey50") +
  # Add text labels for Framingham study
  geom_text(data = subset(clustered_data_unweighted, Study_Population == "Framingham"),
            aes(label = merged_StudyID),
            nudge_y = 0.1, size = 4, color = "blue") +
  # Add text labels for UK-BIOBANK study
  geom_text(data = subset(clustered_data_unweighted, Study_Population == "UK_B"),
            aes(label = merged_StudyID),
            nudge_y = 0.1, size = 4, color = "#C98B9A") +
  # Add text labels for TURKANA_pastoralist ancestry
  geom_text(data = subset(clustered_data_unweighted, Study_Population == "TURKANA_pastoralist"),
            aes(label = merged_StudyID),
            nudge_y = 0.1, size = 4, color = "#8B0000") +
  labs(title = "K-means Clusters Colored by Select Populations",
       x = "Cholesterol Level (scaled)",
       y = "LDL/HDL Ratio (scaled)") +
  # Define color mapping
  scale_color_manual(
    values = c(
      "Circumpolar"         = "#FFB319",
      "TURKANA_pastoralist" = "#8B0000",
      "Framingham"          = "blue",
      "UK_B" = "#C98B9A"
    ),
    na.value = "black"
  ) +
  theme_minimal()
# Save the plot as a PNG file
ggsave(filename = "Clustered_Ancestries_KMEANS.png", plot = plot_unweighted,
       width = 10, height = 8, dpi = 300, bg = "white")

# Read in the data from the CSV file
##run Python ANOVA_AYANA_META.py to get the below file
# Read data
data <- read.csv("/Users/bm0211/RegEx/Water_Energy/anova_results_vs_turkana.csv")
head(data)
dim(data)
colnames(data)
# Arrange by SMD and create unique Population_ID
data <- data %>%
  arrange(Standardized_Mean_Difference) %>%
  group_by(Population) %>%
  mutate(Population_ID = paste(Population, row_number(), sep = "_")) %>%
  ungroup() %>%
  mutate(Population_ID = factor(Population_ID, levels = Population_ID))

# Define effect size categories based on Cohen's d and confidence intervals
data <- data %>%
  mutate(
    Effect_Size_Category = case_when(
      CI_Lower <= 0 & CI_Upper >= 0 ~ "Not Significant", # CI includes 0
      abs(Standardized_Mean_Difference) < 0.2 ~ "Not Significant", # Safety check
      abs(Standardized_Mean_Difference) < 0.5 ~ "Small Effect", # Small effect
      abs(Standardized_Mean_Difference) < 0.8 ~ "Medium Effect", # Medium effect
      TRUE ~ "Large Effect" # Large effect
    )
  )

# Plot 1: Standardized Mean Differences by Population (colored by effect size)
p1 <- ggplot(data, aes(x = Standardized_Mean_Difference, y = Population_ID, color = Effect_Size_Category)) +
  geom_point(size = 1.5) +
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  geom_vline(xintercept = 0.8, linetype = "dashed", color = "#CC5C0E", alpha = 0.5) + # Large effect threshold
  geom_vline(xintercept = -0.8, linetype = "dashed", color = "#CC5C0E", alpha = 0.5) +
  labs(
    title = "Standardized Mean Differences (SMD) of 
    HDL Distribution between Turkana and other Population",
    x = "Standardized Mean Difference (Cohen's d)",
    y = "Population",
    color = "Effect Size"
  ) +
  theme_classic() +   # Removes grid lines and gray background
  scale_color_manual(values = c(
    "Not Significant" = "gray",
    "Small Effect" = "#A2B2D9",
    "Medium Effect" = "#E1AF00",
    "Large Effect" = "#CC5C0E"
  ))

# Save Plot 1 as a PNG file
ggsave("Effect_Size_with_CI_CohensD_ANOVA.png", plot = p1, width = 10, height = 6, dpi = 300)

# Summary of effect sizes (including "Not Significant")
effect_size_summary <- data %>%
  group_by(Effect_Size_Category) %>%
  summarize(Count = n()) %>%
  mutate(Percentage = round(100 * Count / sum(Count), 1))

# Plot 2: Pie chart for effect size categories
p2 <- ggplot(effect_size_summary, aes(x = "", y = Percentage, fill = Effect_Size_Category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  labs(
    title = "Proportion of Effect Size Categories",
    fill = "Effect Size"
  ) +
  geom_text(aes(label = paste0(Percentage, "%")), 
            position = position_stack(vjust = 0.5), 
            color = "white", size = 5) +
  theme_void() +
  scale_fill_manual(values = c(
    "Not Significant" = "gray",
    "Small Effect" = "#A2B2D9",
    "Medium Effect" = "#E1AF00",
    "Large Effect" = "#CC5C0E"
  ))
print(p1)
print(p2)
# Save Plot 2 as a PNG file
ggsave("Effect_Size_with_CI_CohensD_PIE.png", plot = p2, width = 10, height = 6, dpi = 300)


