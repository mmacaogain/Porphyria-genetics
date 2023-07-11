setwd("C:/Users/mmaca/Code/git/Porphyria-genetics/DB_updates/gnomad_1/")

library(tidyverse)
library(dplyr)

Scores<-read.csv("../../Scores.csv")

HMBS <- read.csv("HMBS_gnomAD_v2.1.1_ENSG00000256269_2023_07_04_13_25_11.csv")
CPOX <- read.csv("CPOX_gnomAD_v2.1.1_ENSG00000080819_2023_07_04_14_01_44.csv")
PPOX <- read.csv("PPOX_gnomAD_v2.1.1_ENSG00000143224_2023_07_04_14_02_13.csv")

HMBS$Gene <- "HMBS"
matching_rows <- HMBS$Transcript.Consequence %in% Scores$HGVS_nucleotide
HMBS <- HMBS[!matching_rows, ]
CPOX$Gene <- "CPOX"
matching_rows <- CPOX$Transcript.Consequence %in% Scores$HGVS_nucleotide
CPOX <- CPOX[!matching_rows, ]
PPOX$Gene <- "PPOX"
matching_rows <- PPOX$Transcript.Consequence %in% Scores$HGVS_nucleotide
PPOX <- PPOX[!matching_rows, ]

combined_data <- rbind(HMBS, CPOX, PPOX)
# Filter for rows with column value "missense_variant"
filtered_data <- combined_data %>%
  filter(combined_data$VEP.Annotation == "missense_variant")

gnomadupdate <- data.frame(matrix(nrow = nrow(filtered_data)))
gnomadupdate<-gnomadupdate[,-1]
gnomadupdate$GRCh37<-filtered_data$Position
#gnomadupdate$Chromosome<-filtered_data$Chromosome
gnomadupdate$Gene<-filtered_data$Gene
gnomadupdate$HGVS_nucleotide<-filtered_data$Transcript.Consequence
gnomadupdate$HGVS_protein<-filtered_data$Protein.Consequence
gnomadupdate$ClinVar<-filtered_data$ClinVar.Variation.ID
gnomadupdate$Data_source<-"genomAD_v2.1.1"

# Find common column names
common_cols <- intersect(colnames(gnomadupdate), colnames(Scores))

# Subset data frames to common columns
gnomadupdate_subset <- gnomadupdate[, common_cols]
Scores_subset <- Scores[, common_cols]

# Concatenate the data frames vertically
concatenated_df <- rbind(Scores_subset, gnomadupdate_subset)

# Find the duplicate values in the column
duplicate_values <- concatenated_df$HGVS_protein [duplicated(concatenated_df$HGVS_protein )]

# Print the duplicate values
print(duplicate_values) # will need to look into this in next revision

#Let's rename the original Scores file
file.rename("../../Scores.csv", "../../Scores_old_11072023.csv")

#Merge original with new genomad variants
merged_df_2 <- left_join(concatenated_df, Scores, by = "GRCh37")

#This creates some columns we don't want so remove
merged_df_3 <- select(merged_df_2, -ends_with(".y"))

#but some colnames will still have the '.x' prefix
names(merged_df_3) <- sub("\\.x$", "", names(merged_df_3))
colnames(merged_df_3)

#Write the new updated file to WD. 
write.csv(merged_df_3, "../../Scores.csv", row.names = FALSE)
#Note that the genomic coordinates had to be manually entered in excel. Yep. Excel! also deleted 2x duplicate rows


