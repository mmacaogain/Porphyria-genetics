library(tidyverse)

# Load and filter the data
final_data <- read.csv("gnomad_variant_interpretation/Scores_ANNOVAR_2024_MASTER_WvL_functional.csv") %>%
  filter(Func.refGeneWithVer == "exonic") %>%
  select(-Func.refGeneWithVer, -SJHConsScore, -CONDEL, -HGVS_protein_dual)

# Convert relevant columns to numeric
final_data <- final_data %>%
  mutate(
    AlphaMissense_score = as.numeric(AlphaMissense_score),
    BayesDel_addAF_score_0.07 = as.numeric(BayesDel_addAF_score_0.07),
    REVEL = as.numeric(REVEL),
    SIFT = as.numeric(SIFT),
    MA = as.numeric(MA)
  )

# Define the thresholds
thresholds <- list(
  AlphaMissense_score = 0.5,  # AlphaMissense > 0.5
  BayesDel_addAF_score_0.07 = 0.07,  # BayesDel > 0.07
  REVEL = 0.7,  # REVEL > 0.7
  SIFT = 0.05,  # SIFT < 0.05
  MA = 3.5  # MA > 3.5
)

# Create binary columns to indicate whether each predictor passes the threshold
final_data <- final_data %>%
  mutate(
    AlphaMissense_pass = ifelse(AlphaMissense_score > thresholds$AlphaMissense_score, 1, 0),
    BayesDel_pass = ifelse(BayesDel_addAF_score_0.07 > thresholds$BayesDel_addAF_score_0.07, 1, 0),
    REVEL_pass = ifelse(REVEL > thresholds$REVEL, 1, 0),
    SIFT_pass = ifelse(SIFT <= thresholds$SIFT, 1, 0),  # SIFT is pathogenic if <= threshold
    MA_pass = ifelse(MA > thresholds$MA, 1, 0)
  )

# Calculate the number of valid predictors (non-NA values) for each row
final_data <- final_data %>%
  mutate(
    Valid_Predictors_Count = rowSums(!is.na(select(., AlphaMissense_score, BayesDel_addAF_score_0.07, REVEL, SIFT, MA)))
  )

# Calculate the Consensus Score by summing binary columns and dividing by the number of valid predictors
final_data <- final_data %>%
  mutate(
    Consensus_Score = rowSums(select(., AlphaMissense_pass, BayesDel_pass, REVEL_pass, SIFT_pass, MA_pass), na.rm = TRUE) / Valid_Predictors_Count * 100
  )

# Apply a consensus cutoff (e.g., 60%) to classify variants
consensus_cutoff <- 60  # You can adjust this value
final_data <- final_data %>%
  mutate(Pathogenic_Classification = ifelse(Consensus_Score >= consensus_cutoff, "Pathogenic", "Benign"))

# Print the first few rows to verify
print(head(final_data))

# Filter data to include only the HMBS gene for further analysis
roc_analysis_hmbs <- final_data %>%
  filter(Gene == "HMBS")

# Print the first few rows of roc_analysis_hmbs to verify
print(head(roc_analysis_hmbs))

# List of predictors to evaluate (only the specified ones)
predictors <- c('AlphaMissense_score', 'BayesDel_addAF_score_0.07', 'REVEL', 'SIFT', 'MA', 'Consensus_Score')

# # Ensure the necessary predictors are numeric
# roc_analysis_hmbs$AlphaMissense_score <- as.numeric(roc_analysis_hmbs$AlphaMissense_score)
# roc_analysis_hmbs$BayesDel_addAF_score_0.07 <- as.numeric(roc_analysis_hmbs$BayesDel_addAF_score_0.07)
# roc_analysis_hmbs$REVEL <- as.numeric(roc_analysis_hmbs$REVEL)
# roc_analysis_hmbs$SIFT <- as.numeric(roc_analysis_hmbs$SIFT)
# roc_analysis_hmbs$MA <- as.numeric(roc_analysis_hmbs$MA)
# roc_analysis_hmbs$Consensus_Score <- as.numeric(roc_analysis_hmbs$Consensus_Score)

# Ensure functional data is correctly classified as 0 (benign) and 1 (deleterious)
# Assuming the 'Functional' column is used as ground truth:
# "WT" (wild type) is considered 0 (benign), and "impaired" is considered 1 (deleterious)
roc_analysis_hmbs$Functional <- ifelse(roc_analysis_hmbs$in_vitro_HMBS_activity == "WT", 0, ifelse(roc_analysis_hmbs$in_vitro_HMBS_activity == "Impaired", 1, NA))

# Remove rows with NA in the 'Functional' column
roc_analysis_hmbs <- roc_analysis_hmbs[!is.na(roc_analysis_hmbs$Functional), ]

# Colors for each predictor
colors <- c('AlphaMissense_score' = "blue", 
            'BayesDel_addAF_score_0.07' = "green", 
            'REVEL' = "red", 
            'SIFT' = "purple", 
            'MA' = "orange",
            'Consensus_Score' = 'black')

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

