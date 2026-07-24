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
#file.rename("../../Scores.csv", "../../Scores_old_11072023.csv")

#Merge original with new genomad variants
merged_df_2 <- left_join(concatenated_df, Scores, by = "GRCh37")

#This creates some columns we don't want so remove
merged_df_3 <- select(merged_df_2, -ends_with(".y"))

#but some colnames will still have the '.x' prefix
names(merged_df_3) <- sub("\\.x$", "", names(merged_df_3))
colnames(merged_df_3)

#Write the new updated file to WD. 
#write.csv(merged_df_3, "../../Scores.csv", row.names = FALSE)
#Note that the genomic coordinates had to be manually entered in excel. Yep. Excel! also deleted 2x duplicate rows

#FECH UROD####

# Read in old Scores file (contains only PPOX, CPOX, HMBS)
Scores <- read.csv("../../Scores.csv", stringsAsFactors = FALSE)

# Read new gnomAD variant files for FECH and UROD
FECH <- read.csv("FECH_gnomAD_v2.1.1_ENSG00000066926_2025_03_04_16_55_35.csv", stringsAsFactors = FALSE) %>%
  select(-gnomAD.ID, everything()) # Drop `gnomAD ID` if it exists

UROD <- read.csv("UROD_gnomAD_v2.1.1_ENSG00000126088_2025_03_04_16_58_37.csv", stringsAsFactors = FALSE) %>%
  select(-gnomAD.ID, everything())

# Add Gene column to distinguish datasets
FECH$Gene <- "FECH"
UROD$Gene <- "UROD"

# Filter out variants that are already in Scores (prevents duplicates)
#matching_FECH <- FECH$Transcript.Consequence %in% Scores$HGVS_nucleotide
#FECH <- FECH[!matching_FECH, ]
#matching_UROD <- UROD$Transcript.Consequence %in% Scores$HGVS_nucleotide
#UROD <- UROD[!matching_UROD, ]

# Combine FECH and UROD variants
combined_data <- rbind(FECH, UROD)

# Filter only missense variants
filtered_data <- combined_data %>%
  filter(VEP.Annotation == "missense_variant")

# Create a new data frame for gnomadupdate
gnomadupdate <- data.frame(matrix(nrow = nrow(filtered_data)))
gnomadupdate <- gnomadupdate[,-1] # Remove the empty column

# Select relevant columns
gnomadupdate$GRCh37 <- filtered_data$Position
gnomadupdate$Gene <- filtered_data$Gene
gnomadupdate$HGVS_nucleotide <- filtered_data$Transcript.Consequence
gnomadupdate$HGVS_protein <- filtered_data$Protein.Consequence
gnomadupdate$ClinVarID <- filtered_data$ClinVar.Variation.ID
gnomadupdate$Variant_source <- "gnomAD_v2.1.1"

# Ensure column consistency with Scores before merging
common_cols <- intersect(colnames(gnomadupdate), colnames(Scores))

# Subset data frames to common columns
gnomadupdate_subset <- gnomadupdate[, common_cols]
Scores_subset <- Scores[, common_cols]

# Append new variants to Scores
concatenated_df <- rbind(Scores_subset, gnomadupdate_subset)

# Find duplicate values in the HGVS_protein column
duplicate_values <- concatenated_df$HGVS_protein[duplicated(concatenated_df$HGVS_protein)]

# Print duplicate values for review
print(duplicate_values)

#Rename the old Scores file before overwriting
file.rename("../../Scores.csv", "../../Scores_old_13032025.csv")

#Merge old and new variants by GRCh37 position
merged_df_2 <- left_join(concatenated_df, Scores, by = "GRCh37")

#Remove unwanted columns ending in ".y"
merged_df_3 <- select(merged_df_2, -ends_with(".y"))

#Rename columns by removing ".x" suffixes
names(merged_df_3) <- sub("\\.x$", "", names(merged_df_3))

#Note that the genomic will now have to be manually entered in EXCEL! also deleted 2x duplicate rows

# Write updated Scores file
write.csv(merged_df_3, "../../Scores_FECH_UROD.csv", row.names = FALSE)

