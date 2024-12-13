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
  select(-Func.refGeneWithVer, -SJHConsScore, -in_vitro_HMBS_activity, 
         -CONDEL, -SJHConsScore, -HGVS_protein_dual, -DiscovEHR
  )

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
    Consensus_Score = ceiling(rowSums(select(., PP2_pass, SIFT_pass, Consurf_pass, MA_pass, PROVEAN_pass, FATHMM_pass, REVEL_pass, CADD_pass, AlphaMissense_pass, BayesDel_pass), na.rm = TRUE) / Valid_Predictors_Count * 100)
  )

#Apply a cut off of at >3 predictor scores and drop unneccessary columns
final_data_cns <- final_data %>%
  filter(Valid_Predictors_Count > 3)%>%
  select(-Valid_Predictors_Count, -SIFT_pass, -PP2_pass, -Consurf_pass, -MA_pass, -PROVEAN_pass, -FATHMM_pass, -REVEL_pass, -CADD_pass, -AlphaMissense_pass, -BayesDel_pass)

#View the first few rows of the final dataset
head(final_data_cns)

final_data_cns <- final_data_cns %>% rename(BayesDel = BayesDel_addAF_score_0.07)
final_data_cns <- final_data_cns %>% rename(AlphaMissense = AlphaMissense_score)


# Optionally, save the filtered data to a new CSV file
write.csv(final_data_cns, "../../Scores.csv", row.names = FALSE)


