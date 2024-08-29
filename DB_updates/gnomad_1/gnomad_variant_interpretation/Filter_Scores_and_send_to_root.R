setwd("C:/Users/mmaca/Code/git/Porphyria-genetics/DB_updates/gnomad_1/")

library(tidyverse)
library(dplyr)

# Read in the CSV file, filter rows, and drop specified columns in one block
final_data <- read.csv("gnomad_variant_interpretation/Scores_ANNOVAR_2024_MASTER.csv") %>%
  filter(Func.refGeneWithVer == "exonic") %>%
  select(-Func.refGeneWithVer, -SJHConsScore, -in_vitro_HMBS_activity, -SJHConsScore, -HGVS_protein_dual
)

# View the first few rows of the final dataset
head(final_data)

# Optionally, save the filtered data to a new CSV file
write.csv(final_data, "../../Scores.csv", row.names = FALSE)
