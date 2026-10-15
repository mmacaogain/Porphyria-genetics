#Logistic regression####
library(pROC)
library(dplyr)

# Step 1: Define default and adjusted cut-offs
default_cutoffs <- list(
  AlphaMissense_score = 0.5,     # Default cut-off for AlphaMissense
  BayesDel_addAF_score_0.07 = 0.07, # Default cut-off for BayesDel
  REVEL = 0.7,                   # Default cut-off for REVEL
  SIFT = 0.05,                   # Default cut-off for SIFT
  MA = 3.5,                      # Default cut-off for MA
  Consensus_Score = 50            # Default cut-off for Consensus_Score (example)
)

adjusted_cutoffs <- list(
  AlphaMissense_score = 0.76,     # Adjusted cut-off for AlphaMissense
  BayesDel_addAF_score_0.07 = 0.44, # Adjusted cut-off for BayesDel
  REVEL = 0.9,                    # Adjusted cut-off for REVEL
  SIFT = 0.03,                    # Adjusted cut-off for SIFT
  MA = 2.54,                      # Adjusted cut-off for MA
  Consensus_Score = 50            # Adjusted cut-off for Consensus_Score (example)
)

# Step 2: Apply the default cut-offs to create binary variables for each predictor
roc_analysis_hmbs <- roc_analysis_hmbs %>%
  mutate(
    AlphaMissense_pass_default = ifelse(AlphaMissense_score > default_cutoffs$AlphaMissense_score, 1, 0),
    BayesDel_pass_default = ifelse(BayesDel_addAF_score_0.07 > default_cutoffs$BayesDel_addAF_score_0.07, 1, 0),
    REVEL_pass_default = ifelse(REVEL > default_cutoffs$REVEL, 1, 0),
    SIFT_pass_default = ifelse(SIFT <= default_cutoffs$SIFT, 1, 0),  # For SIFT, lower is pathogenic
    MA_pass_default = ifelse(MA > default_cutoffs$MA, 1, 0)
  )

# Step 3: Apply the adjusted cut-offs to create binary variables for each predictor
roc_analysis_hmbs <- roc_analysis_hmbs %>%
  mutate(
    AlphaMissense_pass_adjusted = ifelse(AlphaMissense_score > adjusted_cutoffs$AlphaMissense_score, 1, 0),
    BayesDel_pass_adjusted = ifelse(BayesDel_addAF_score_0.07 > adjusted_cutoffs$BayesDel_addAF_score_0.07, 1, 0),
    REVEL_pass_adjusted = ifelse(REVEL > adjusted_cutoffs$REVEL, 1, 0),
    SIFT_pass_adjusted = ifelse(SIFT <= adjusted_cutoffs$SIFT, 1, 0),
    MA_pass_adjusted = ifelse(MA > adjusted_cutoffs$MA, 1, 0)
  )

# Step 4: Train logistic regression models for both default and adjusted cut-offs

# Logistic Regression with Default Cut-offs
logistic_model_default <- glm(Functional ~ AlphaMissense_pass_default + BayesDel_pass_default + 
                                REVEL_pass_default + SIFT_pass_default + MA_pass_default, 
                              data = roc_analysis_hmbs, family = binomial)

# Logistic Regression with Adjusted Cut-offs
logistic_model_adjusted <- glm(Functional ~ AlphaMissense_pass_adjusted + BayesDel_pass_adjusted + 
                                 REVEL_pass_adjusted + SIFT_pass_adjusted + MA_pass_adjusted, 
                               data = roc_analysis_hmbs, family = binomial)

# Step 5: Predict probabilities for both models

# Logistic Regression Predictions - Default
logistic_probs_default <- predict(logistic_model_default, roc_analysis_hmbs, type = "response")

# Logistic Regression Predictions - Adjusted
logistic_probs_adjusted <- predict(logistic_model_adjusted, roc_analysis_hmbs, type = "response")

# Step 6: Calculate ROC Curves for both models
logistic_roc_default <- roc(roc_analysis_hmbs$Functional, logistic_probs_default)
logistic_roc_adjusted <- roc(roc_analysis_hmbs$Functional, logistic_probs_adjusted)

# Step 7: Print AUC values for both models
cat("Logistic Regression Default AUC:", auc(logistic_roc_default), "\n")
cat("Logistic Regression Adjusted AUC:", auc(logistic_roc_adjusted), "\n")

# Step 8: Plot both ROC curves on the same plot (Default in Blue, Adjusted in Red)
plot(logistic_roc_default, col = "blue", lwd = 2, main = "ROC Curve: Default vs Adjusted Cut-offs")
plot(logistic_roc_adjusted, col = "red", lwd = 2, add = TRUE)  # Add Adjusted ROC to the same plot
legend("bottomright", legend = c("Default Cut-offs", "Adjusted Cut-offs"), col = c("blue", "red"), lwd = 2)

# Perform DeLong's test to compare the two ROC curves
roc_test <- roc.test(logistic_roc_default, logistic_roc_adjusted)

# Print the results of the test
print(roc_test)

#Random Forrest####
library(randomForest)
library(pROC)

# Step 1: Train Random Forest Model with Default Cut-offs
rf_model_default <- randomForest(as.factor(Functional) ~ AlphaMissense_pass_default + BayesDel_pass_default + 
                                   REVEL_pass_default + SIFT_pass_default + MA_pass_default, 
                                 data = roc_analysis_hmbs, ntree = 500)

# Step 2: Train Random Forest Model with Adjusted Cut-offs
rf_model_adjusted <- randomForest(as.factor(Functional) ~ AlphaMissense_pass_adjusted + BayesDel_pass_adjusted + 
                                    REVEL_pass_adjusted + SIFT_pass_adjusted + MA_pass_adjusted, 
                                  data = roc_analysis_hmbs, ntree = 500)

# Step 3: Predict Probabilities for both models
rf_probs_default <- predict(rf_model_default, roc_analysis_hmbs, type = "prob")[, 2]  # Probability for class 1 (Pathogenic)
rf_probs_adjusted <- predict(rf_model_adjusted, roc_analysis_hmbs, type = "prob")[, 2]  # Probability for class 1 (Pathogenic)

# Step 4: Calculate ROC Curves for both models
rf_roc_default <- roc(roc_analysis_hmbs$Functional, rf_probs_default)
rf_roc_adjusted <- roc(roc_analysis_hmbs$Functional, rf_probs_adjusted)

# Step 5: Print AUC values for both models
cat("Random Forest Default AUC:", auc(rf_roc_default), "\n")
cat("Random Forest Adjusted AUC:", auc(rf_roc_adjusted), "\n")

# Step 6: Perform DeLong's test to compare the two ROC curves
rf_roc_test <- roc.test(rf_roc_default, rf_roc_adjusted)

# Print the results of the test
print(rf_roc_test)

# Step 7: Plot both ROC curves on the same plot (Default in Blue, Adjusted in Red)
plot(rf_roc_default, col = "blue", lwd = 2, main = "Random Forest ROC: Default vs Adjusted Cut-offs")
plot(rf_roc_adjusted, col = "red", lwd = 2, add = TRUE)  # Add Adjusted ROC to the same plot
legend("bottomright", legend = c("Default Cut-offs", "Adjusted Cut-offs"), col = c("blue", "red"), lwd = 2)
