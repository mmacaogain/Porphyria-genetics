# List of predictors to evaluate (only the specified ones)
predictors <- c('AlphaMissense_score', 'BayesDel_addAF_score_0.07', 'REVEL', 'SIFT', 'MA', 'Consensus_Score')

# Ensure the necessary predictors are numeric
roc_analysis_hmbs$AlphaMissense_score <- as.numeric(roc_analysis_hmbs$AlphaMissense_score)
roc_analysis_hmbs$BayesDel_addAF_score_0.07 <- as.numeric(roc_analysis_hmbs$BayesDel_addAF_score_0.07)
roc_analysis_hmbs$REVEL <- as.numeric(roc_analysis_hmbs$REVEL)
roc_analysis_hmbs$SIFT <- as.numeric(roc_analysis_hmbs$SIFT)
roc_analysis_hmbs$MA <- as.numeric(roc_analysis_hmbs$MA)
roc_analysis_hmbs$Consensus_Score <- as.numeric(roc_analysis_hmbs$Consensus_Score)

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
