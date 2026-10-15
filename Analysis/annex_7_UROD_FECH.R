setwd("C:/Users/mmaca/Code/git/Porphyria-genetics/Analysis/")
library(reshape2)
library(knitr)
library(randomForest)
library(caret)
library(pROC)
library(dplyr)
library(tidyverse)

# Load and filter the data
final_data <- read.csv("../DB_updates/gnomad_1/gnomad_variant_interpretation/Scores_ANNOVAR_2025_MASTER_WvL_functional_UROD_FECH.csv") %>%
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
  Median_Score_PPOX = c(0.220000, 0.142613, 1.952500, 0.634000, 0.091000),
  Median_Score_UROD = c(0.167000, 0.177395, 2.195000, 0.680000, 0.049000),
  Median_Score_FECH = c(0.2730000, 0.2698835, 1.9450000, 0.6975000, 0.0405000)
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
  PPOX_Cutoff = HMBS_Cutoff * (summary_stats$Median_Score_PPOX / summary_stats$Median_Score_HMBS),
  
  # Adjust UROD cutoffs based on the ratio of UROD median to HMBS median
  UROD_Cutoff = HMBS_Cutoff * (summary_stats$Median_Score_UROD / summary_stats$Median_Score_HMBS),
  
  # Adjust FECH cutoffs based on the ratio of FECH median to HMBS median
  FECH_Cutoff = HMBS_Cutoff * (summary_stats$Median_Score_FECH / summary_stats$Median_Score_HMBS),
  
)


# Print the adjusted cutoffs
print(adjusted_cutoffs)

# Step 3: Apply the cutoffs to classify variants based on the adjusted cutoffs
long_data <- long_data %>%
  left_join(adjusted_cutoffs, by = "Predictor") %>%
  mutate(Cutoff = case_when(
    Gene == "HMBS" ~ HMBS_Cutoff,
    Gene == "CPOX" ~ CPOX_Cutoff,
    Gene == "PPOX" ~ PPOX_Cutoff,
    Gene == "UROD" ~ UROD_Cutoff,
    Gene == "FECH" ~ FECH_Cutoff
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
  mutate(Gene = factor(Gene, levels = c("HMBS", "PPOX", "CPOX", "UROD", "FECH")))

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
