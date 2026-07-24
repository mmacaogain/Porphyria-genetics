setwd("C:/Users/mmaca/Code/git/Porphyria-genetics/DB_updates/gnomad_1/")

library(tidyverse)
library(dplyr)
library(pROC)

#Read in the CSV file, filter rows, and drop specified columns in one block
# final_data <- read.csv("gnomad_variant_interpretation/Scores_ANNOVAR_2024_MASTER.csv") %>%
#   filter(Func.refGeneWithVer == "exonic") %>%
#   select(-Func.refGeneWithVer, -SJHConsScore, #-in_vitro_HMBS_activity, 
#          -CONDEL, -SJHConsScore, -HGVS_protein_dual
#   )

final_data <- read.csv("gnomad_variant_interpretation/Scores_ANNOVAR_2024_MASTER_WvL_functional.csv") %>%
  filter(Func.refGeneWithVer == "exonic") %>%
  select(-Func.refGeneWithVer, -SJHConsScore, #-in_vitro_HMBS_activity, 
         -CONDEL, -SJHConsScore, -HGVS_protein_dual, -DiscovEHR
  )
#FOR UPDATES ADD THIS "gnomad_variant_interpretation/Scores_ANNOVAR_2025_MASTER_WvL_functional_UROD_FECH.csv"

# Define the thresholds for each predictor
thresholds <- c(
  PP2 = 0.908,  # Threshold for PolyPhen-2
  SIFT = 0.05,  # Threshold for SIFT
  Consurf = 7,  # Threshold for Consurf
  MA = 3.5,  # Threshold for MutationAssessor (MA)
  PROVEAN = -2.5,  # Threshold for PROVEAN
  FATHMM = -2.0,  # Threshold for FATHMM
  REVEL = 0.7,  # Threshold for REVEL
  CADD = 25,  # Threshold for CADD (PHRED score)
  AlphaMissense_score = 0.5,  # Threshold for AlphaMissense
  #DiscovEHR = 1,  # DiscovEHR is treated as 1 if not NA, 0 if NA
  BayesDel_addAF_score_0.07 = 0.07  # Threshold for BayesDel_addAF
)

# Create binary columns to indicate whether each predictor passes the threshold
final_data <- final_data %>%
  mutate(
    PP2_pass = ifelse(PP2 > thresholds["PP2"], 1, 0),  # PolyPhen-2
    SIFT_pass = ifelse(SIFT <= thresholds["SIFT"], 1, 0),  # SIFT is pathogenic if <= threshold
    Consurf_pass = ifelse(Consurf >= thresholds["Consurf"], 1, 0),  # Consurf
    MA_pass = ifelse(MA > thresholds["MA"], 1, 0),  # MutationAssessor
    PROVEAN_pass = ifelse(PROVEAN <= thresholds["PROVEAN"], 1, 0),  # PROVEAN
    FATHMM_pass = ifelse(FATHMM <= thresholds["FATHMM"], 1, 0),  # FATHMM
    REVEL_pass = ifelse(REVEL > thresholds["REVEL"], 1, 0),  # REVEL
    CADD_pass = ifelse(CADD >= thresholds["CADD"], 1, 0),  # CADD
    AlphaMissense_pass = ifelse(AlphaMissense_score > thresholds["AlphaMissense_score"], 1, 0),  # AlphaMissense
    #DiscovEHR_pass = ifelse(!is.na(DiscovEHR) & DiscovEHR == thresholds["DiscovEHR"], 1, 0),  # DiscovEHR
    BayesDel_pass = ifelse(BayesDel_addAF_score_0.07 > thresholds["BayesDel_addAF_score_0.07"], 1, 0)  # BayesDel
  )

# Calculate the number of valid predictors (non-NA values) for each row
final_data <- final_data %>%
  mutate(
    Valid_Predictors_Count = rowSums(!is.na(select(., PP2, SIFT, Consurf, MA, PROVEAN, FATHMM, REVEL, CADD, AlphaMissense_score, BayesDel_addAF_score_0.07)))
  )

# Calculate the Consensus Score by summing binary columns and dividing by the number of valid predictors
final_data <- final_data %>%
  mutate(
    Consensus_Score = rowSums(select(., PP2_pass, SIFT_pass, Consurf_pass, MA_pass, PROVEAN_pass, FATHMM_pass, REVEL_pass, CADD_pass, AlphaMissense_pass, BayesDel_pass), na.rm = TRUE) / Valid_Predictors_Count * 100
  )

# Apply a consensus cutoff (e.g., 60%) to classify variants
consensus_cutoff <- 60  # You can adjust this value
final_data <- final_data %>%
  mutate(Pathogenic_Classification = ifelse(Consensus_Score >= consensus_cutoff, "Pathogenic", "Benign"))

# Print the first few rows to verify
print(head(final_data))

final_data_cns <- final_data %>%
   filter(Valid_Predictors_Count > 3)%>%
  select(-Valid_Predictors_Count
  )

# #View the first few rows of the final dataset
# head(final_data_cns)
# 
# # Optionally, save the filtered data to a new CSV file
# write.csv(final_data_cns, "../../Scores.csv", row.names = FALSE)

roc_analysis <- final_data_cns

# Load the data (assuming you have already loaded your dataset)
# Filter the data to include only HMBS gene
roc_analysis_hmbs <- roc_analysis[roc_analysis$Gene == "HMBS", ]

# Ensure functional data is correctly classified as 0 (benign) and 1 (deleterious)
# Assuming the 'Functional' column is used as ground truth:
# "WT" (wild type) is considered 0 (benign), and "impaired" is considered 1 (deleterious)
roc_analysis_hmbs$Functional <- ifelse(roc_analysis_hmbs$in_vitro_HMBS_activity == "WT", 0, ifelse(roc_analysis_hmbs$in_vitro_HMBS_activity == "Impaired", 1, NA))

# Remove rows with NA in the 'Functional' column
roc_analysis_hmbs <- roc_analysis_hmbs[!is.na(roc_analysis_hmbs$Functional), ]

# List of predictors to evaluate (only the specified ones)
predictors <- c('PP2','SIFT', 'Consurf',
                'MA',
                'PROVEAN', 
                'FATHMM',  
                'REVEL', 
                'CADD', 
                'AlphaMissense_score', 
                'BayesDel_addAF_score_0.07', 
                'Consensus_Score')

# Colors for each predictor
colors <- c(
  'PP2' = "darkblue", 
  'SIFT' = "purple", 
  'Consurf' = "lightblue", 
  'MA' = "orange", 
  'PROVEAN' = "pink", 
  'FATHMM' = "brown", 
  'REVEL' = "red", 
  'CADD' = "cyan", 
  'AlphaMissense_score' = "blue", 
  'DiscovEHR' = "darkgreen", 
  'BayesDel_addAF_score_0.07' = "green", 
  'Consensus_Score' = 'black'
)

# Initialize a list to store AUC values
auc_values <- list()

# Create an empty plot for ROC curves
plot(NULL, xlim = c(0, 1), ylim = c(0, 1), xlab = "False Positive Rate", ylab = "True Positive Rate",
     main = "ROC Curves for Selected Predictors")

# Add diagonal line for reference
abline(a = 0, b = 1, lty = 2, col = "gray")

# Generate and plot ROC curve for each predictor
for (predictor in predictors) {
  roc_obj <- roc(roc_analysis_hmbs$Functional, roc_analysis_hmbs[[predictor]])
  auc_values[[predictor]] <- auc(roc_obj)
  
  # Plot sensitivity (y) vs. specificity (1 - x)
  lines(1 - roc_obj$specificities, roc_obj$sensitivities, col = colors[predictor], lwd = 2)
}

# Add a legend to the plot
legend("bottomright", legend = predictors, col = colors[predictors], lwd = 2)

# Display AUC values
auc_values

# If you want to save the plots, you can wrap the plotting code with pdf() and dev.off()
# pdf("roc_curves.pdf")
# for (predictor in predictors) {
#   generate_roc(predictor, hmbs_data)
# }
# dev.off()

library(pROC)

# Assuming hmbs_data_clean is the cleaned data with all predictors as numeric

# Initialize a list to store optimal cutoffs
optimal_cutoffs <- list()

# Generate ROC curve and determine the optimal cutoff for each predictor
for (predictor in predictors) {
  roc_obj <- roc(roc_analysis_hmbs$Functional, roc_analysis_hmbs[[predictor]])
  
  # Calculate Youden's J statistic for each point on the ROC curve
  youden_index <- roc_obj$sensitivities + roc_obj$specificities - 1
  
  # Find the index of the maximum Youden's J statistic
  optimal_index <- which.max(youden_index)
  
  # Get the corresponding cutoff
  optimal_cutoff <- roc_obj$thresholds[optimal_index]
  
  # Store the optimal cutoff
  optimal_cutoffs[[predictor]] <- optimal_cutoff
  
  # Plot ROC curve
  plot(roc_obj, main = paste("ROC Curve for", predictor))
  abline(h = roc_obj$sensitivities[optimal_index], col = "blue", lty = 2)
  abline(v = roc_obj$specificities[optimal_index], col = "blue", lty = 2)
  text(0.5, 0.2, paste("Optimal cutoff:", round(optimal_cutoff, 3)), col = "red")
}

# Display the optimal cutoffs for each predictor
optimal_cutoffs

# List of predictors to evaluate (only the specified ones)
predictors <- c('AlphaMissense_score', 'BayesDel_addAF_score_0.07', 'REVEL', 'SIFT', 'MA')

# Ensure the necessary predictors are numeric
roc_analysis_hmbs$AlphaMissense_score <- as.numeric(roc_analysis_hmbs$AlphaMissense_score)
roc_analysis_hmbs$BayesDel_addAF_score_0.07 <- as.numeric(roc_analysis_hmbs$BayesDel_addAF_score_0.07)
roc_analysis_hmbs$REVEL <- as.numeric(roc_analysis_hmbs$REVEL)
roc_analysis_hmbs$SIFT <- as.numeric(roc_analysis_hmbs$SIFT)
roc_analysis_hmbs$MA <- as.numeric(roc_analysis_hmbs$MA)

# Function to generate ROC curve for each predictor
generate_roc <- function(predictor, data) {
  roc_obj <- roc(data$Functional, data[[predictor]])
  plot(roc_obj, main = paste("ROC Curve for", predictor), col = "blue")
  auc <- auc(roc_obj)
  print(paste("AUC for", predictor, ":", auc))
  return(auc)
}

# Initialize a list to store AUC values
auc_values <- list()

# Generate and plot ROC curve for each predictor
par(mfrow = c(2, 3))  # Adjust the layout to fit all plots
for (predictor in predictors) {
  auc_values[[predictor]] <- generate_roc(predictor, roc_analysis_hmbs)
}

# Reset layout
par(mfrow = c(1, 1))

# Display AUC values
auc_values

