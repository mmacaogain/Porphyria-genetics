# Assuming DeLong's test has already been run
# DeLong's test results (replace these with your actual values)
delong_p_value <- roc_test$p.value
delong_ci <- roc_test$ci  # Confidence interval from DeLong's test
auc_default <- auc(logistic_roc_default)  # AUC for default model
auc_adjusted <- auc(logistic_roc_adjusted)  # AUC for adjusted model

# Plot the ROC curves
plot(logistic_roc_default, col = "blue", lwd = 2, main = "ROC Curve: Default vs Adjusted Cut-offs")
plot(logistic_roc_adjusted, col = "red", lwd = 2, add = TRUE)
legend("bottomright", legend = c(paste("Default (AUC =", round(auc_default, 3), ")"),
                                 paste("Adjusted (AUC =", round(auc_adjusted, 3), ")")),
       col = c("blue", "red"), lwd = 2)

# Add text to display the DeLong test results
text(0.5, 0.1, labels = paste("DeLong's p-value: ", round(delong_p_value, 4)), adj = 0)
text(0.5, 0.05, labels = paste("95% CI: [", round(delong_ci[1], 3), ", ", round(delong_ci[2], 3), "]"), adj = 0)
