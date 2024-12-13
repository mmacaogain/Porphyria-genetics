setwd("C:/Users/mmaca/Code/git/Porphyria-genetics/DB_updates/gnomad_1/")

library(tidyverse)
library(dplyr)
library(pROC)

# Load and filter the data
final_data <- read.csv("gnomad_variant_interpretation/Scores_ANNOVAR_2024_MASTER_WvL_functional.csv") %>%
  filter(Func.refGeneWithVer == "exonic") %>%
  select(-Func.refGeneWithVer, -SJHConsScore, #-in_vitro_HMBS_activity, 
         -CONDEL, -SJHConsScore, -HGVS_protein_dual, -DiscovEHR
  )

#Convert relevant columns to numeric
final_data <- final_data %>%
  mutate(
    AlphaMissense_score = as.numeric(AlphaMissense_score),
    BayesDel_addAF_score_0.07 = as.numeric(BayesDel_addAF_score_0.07),
    REVEL = as.numeric(REVEL),
    SIFT = as.numeric(SIFT),
    MA = as.numeric(MA)
  )

# Define the thresholds []
thresholds <- list(
  AlphaMissense_score = 0.76,  # AlphaMissense > 0.5 / 0.76
  BayesDel_addAF_score_0.07 = 0.44,  # BayesDel > 0.07 / 0.44
  REVEL = 0.9,  # REVEL > 0.7 / 0.9
  SIFT = 0.03,  # SIFT < 0.05 / 0.03
  MA = 2.54  # MA > 3.5 / 2.54
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

#Initialize a list to store AUC values
auc_values <- list()

#Create an empty plot for ROC curves
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
#legend("bottomright", legend = predictors, col = colors[predictors], lwd = 2)

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

library(tidyverse)

# Step 1: Reshape the data and calculate mean and median for each gene and predictor
long_data <- final_data %>%
  pivot_longer(cols = c(AlphaMissense_score, BayesDel_addAF_score_0.07, REVEL, SIFT, MA),
               names_to = "Predictor",
               values_to = "Score")

# Calculate summary statistics for HMBS, CPOX, and PPOX
summary_stats <- final_data %>%
  pivot_longer(cols = c(AlphaMissense_score, BayesDel_addAF_score_0.07, REVEL, SIFT, MA),
               names_to = "Predictor",
               values_to = "Score") %>%
  group_by(Gene, Predictor) %>%
  summarize(Mean_Score = mean(Score, na.rm = TRUE),
            Median_Score = median(Score, na.rm = TRUE)) %>%
  pivot_wider(names_from = Gene, values_from = c(Mean_Score, Median_Score))

####
# Function to get the percentile rank of a cutoff in the HMBS distribution
get_percentile_rank <- function(cutoff, predictor_column, data) {
  ecdf_func <- ecdf(data[[predictor_column]])  # Empirical CDF function
  percentile_rank <- ecdf_func(cutoff)
  return(percentile_rank)
}

# Function to calculate the adjusted cutoff for a different gene based on percentile rank
adjust_cutoff_for_gene <- function(percentile_rank, predictor_column, data) {
  quantile_func <- quantile(data[[predictor_column]], probs = percentile_rank)
  return(quantile_func)
}

# List of predictor columns
predictor_columns <- c("AlphaMissense_score", "BayesDel_addAF_score_0.07", "REVEL", "SIFT", "MA")

# Create the summary stats tibble with the data provided
summary_stats <- tibble(
  Predictor = c("AlphaMissense_score", "BayesDel_addAF_score_0.07", "MA", "REVEL", "SIFT"),
  Median_Score_CPOX = c(0.264500, 0.280913, 2.307500, 0.721000, 0.032000),
  Median_Score_HMBS = c(0.639500, 0.393788, 2.595000, 0.853000, 0.013000),
  Median_Score_PPOX = c(0.220000, 0.142613, 1.952500, 0.634000, 0.091000)
)

# The original HMBS cutoff values
HMBS_Cutoff <- c(0.761, 0.442201, 2.5425, 0.901, 0.03)

# Adjusted cutoffs for CPOX and PPOX based on the ratio of median scores
adjusted_cutoffs <- tibble(
  Predictor = c("AlphaMissense_score", "BayesDel_addAF_score_0.07", "MA", "REVEL", "SIFT"),
  HMBS_Cutoff = HMBS_Cutoff,
  
  # Adjust CPOX cutoffs based on the ratio of CPOX median to HMBS median
  CPOX_Cutoff = HMBS_Cutoff * (summary_stats$Median_Score_CPOX / summary_stats$Median_Score_HMBS),
  
  # Adjust PPOX cutoffs based on the ratio of PPOX median to HMBS median
  PPOX_Cutoff = HMBS_Cutoff * (summary_stats$Median_Score_PPOX / summary_stats$Median_Score_HMBS)
)

# Print the adjusted cutoffs
print(adjusted_cutoffs)

# Step 3: Apply the cutoffs to classify variants based on the adjusted cutoffs
long_data <- long_data %>%
  left_join(adjusted_cutoffs, by = "Predictor") %>%
  mutate(Cutoff = case_when(
    Gene == "HMBS" ~ HMBS_Cutoff,
    Gene == "CPOX" ~ CPOX_Cutoff,
    Gene == "PPOX" ~ PPOX_Cutoff
  )) %>%
  mutate(Pathogenic_Classification = ifelse(Score >= Cutoff, "Pathogenic", "Benign"))

# Step 4: Visualize the results
ggplot(long_data, aes(x = Gene, y = Score, color = in_vitro_HMBS_activity)) +
  geom_jitter(size = 1, alpha = 0.6, width = .2) +  # Plot points on top of violins
  # Add violin plots behind the points
  geom_boxplot(aes(fill = Gene), alpha = 0.3, color = NA) +  # Adjust alpha for transparency, no outline color
  geom_hline(aes(yintercept = Cutoff), linetype = "dashed", color = "red") +  # Add cutoff lines
  facet_wrap(~ Predictor, scales = "free", nrow = 1) +  # Create a facet for each predictor
  labs(title = "Scores for Each Predictor by Gene with Adjusted Cutoff Lines",
       x = "Gene",
       y = "Score") +
  theme_minimal() +
  theme(strip.text = element_text(size = 12))  # Adjust text size in the facet labels

###alt

medians <- long_data %>%
  group_by(Gene, Predictor) %>%
  summarise(Median_Score = median(Score, na.rm = TRUE)) %>%
  ungroup()

# Step 1: Define the default and optimized cutoffs for each predictor
cutoffs <- tibble(
  Predictor = c("AlphaMissense_score", "BayesDel_addAF_score_0.07", "REVEL", "SIFT", "MA"),
  Default_Cutoff = c(0.5, 0.07, 0.7, 0.05, 3.5),
  Optimized_Cutoff = c(0.76, 0.44, 0.9, 0.03, 2.54)
)

# Step 2: Create the cutoffs data frame specifically for HMBS cutoffs
hmbs_cutoffs <- cutoffs %>%
  rename(Cutoff = Optimized_Cutoff) %>%
  mutate(Gene = "HMBS")

# Step 3: Reorder Gene factor levels (HMBS, PPOX, CPOX)
long_data <- long_data %>%
  mutate(Gene = factor(Gene, levels = c("HMBS", "PPOX", "CPOX")))

# Step 4: Update predictor labels
predictor_labels <- c(
  "AlphaMissense_score" = "Alpha missense",
  "BayesDel_addAF_score_0.07" = "BayesDel",
  "MA" = "MA",
  "REVEL" = "REVEL",
  "SIFT" = "SIFT"
)

# Step 5: Plot the results with violin plots, jittered points, and median bars
ggplot(long_data, aes(x = Gene, y = Score, color = in_vitro_HMBS_activity)) +
  
  # Add violin plots behind the points, removing the legend for the Gene fill
  geom_violin(aes(fill = Gene), alpha = 0.3, color = NA, show.legend = FALSE) +  # No outline color for violins
  
  # Add jittered points
  geom_jitter(size = 1, alpha = 0.6, width = 0.2) +  # Jittered points to show individual data points
  
  # Add dashed red lines for the optimized cutoff values for HMBS
  geom_hline(data = hmbs_cutoffs, aes(yintercept = Cutoff), linetype = "dashed", color = "red") +  # Optimized cutoffs
  
  # Add solid blue lines for the default cutoff values
  geom_hline(data = cutoffs, aes(yintercept = Default_Cutoff), linetype = "solid", color = "blue") +  # Default cutoffs
  
  # Add black bars to represent the medians for each Gene-Predictor combination
  geom_crossbar(data = medians, aes(x = Gene, ymin = Median_Score, ymax = Median_Score, y = Median_Score), 
                color = "black", width = 0.5, fatten = 2) +  # Median bars as horizontal lines
  
  # Facet by predictor with custom labels
  facet_wrap(~ Predictor, scales = "free", nrow = 2, labeller = as_labeller(predictor_labels)) +  # Separate plots for each predictor with new labels
  
  # Labels and theme
  labs(title = "Scores for Each Predictor by Gene with Default and Optimized Cutoff Lines",
       x = "Gene",
       y = "Score") +
  theme_minimal() +
  theme(strip.text = element_text(size = 12))  # Adjust text size in facet labels






######
library(randomForest)
library(pROC)

# Step 2: Train Logistic Regression Model using continuous score columns
logistic_model <- glm(Functional ~ Consensus_Score + BayesDel_addAF_score_0.07 + 
                        REVEL + SIFT + MA + AlphaMissense_score, data = roc_analysis_hmbs, family = binomial)

summary(logistic_model)

# Step 3: Train Random Forest Model using continuous score columns
# Ensure Random Forest knows it's a classification problem by using a binary target
rf_model <- randomForest(as.factor(Functional) ~ AlphaMissense_score + BayesDel_addAF_score_0.07 + 
                           REVEL + SIFT + MA + Consensus_Score, data = roc_analysis_hmbs, ntree = 500)

# Step 4: Predict Probabilities for Logistic Regression and Random Forest
# Logistic Regression Predictions
logistic_probs <- predict(logistic_model, roc_analysis_hmbs, type = "response")

# Random Forest Predictions (class probabilities)
rf_probs <- predict(rf_model, roc_analysis_hmbs, type = "prob")[, 2]  # Probability for class 1 (Pathogenic)

# Step 5: Calculate ROC Curves for Both Models
logistic_roc <- roc(roc_analysis_hmbs$Functional, logistic_probs)
rf_roc <- roc(roc_analysis_hmbs$Functional, rf_probs)

# Step 6: Print AUC values for both models
cat("Logistic Regression AUC:", auc(logistic_roc), "\n")
cat("Random Forest AUC:", auc(rf_roc), "\n")

# Step 7: Plot ROC Curves for Both Models
plot(logistic_roc, col = "blue", lwd = 2, main = "ROC Curve Comparison: Logistic vs Random Forest")
plot(rf_roc, col = "green", lwd = 2, add = TRUE)  # Add random forest ROC to the same plot
legend("bottomright", legend = c("Logistic Regression", "Random Forest"), col = c("blue", "green"), lwd = 2)

library(caret)
library(randomForest)

# Set up cross-validation
set.seed(123)
train_control <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = twoClassSummary)

# Convert Functional to a factor with valid R variable names
roc_analysis_hmbs$Functional <- as.factor(roc_analysis_hmbs$Functional)

# Rename the levels of Functional to something meaningful and valid
levels(roc_analysis_hmbs$Functional) <- c("Benign", "Pathogenic")

# Alternatively, if using Pathogenic_Classification_Binary, apply the same fix
roc_analysis_hmbs$Pathogenic_Classification_Binary <- as.factor(roc_analysis_hmbs$Pathogenic_Classification_Binary)
levels(roc_analysis_hmbs$Pathogenic_Classification_Binary) <- c("Benign", "Pathogenic")


# Train the Random Forest model with cross-validation
rf_cv_model <- train(
  Functional ~ AlphaMissense_score + BayesDel_addAF_score_0.07 + REVEL + SIFT + MA,
  data = roc_analysis_hmbs,
  method = "rf",
  trControl = train_control,
  metric = "ROC"
)

# Print the cross-validated AUC
print(rf_cv_model)

