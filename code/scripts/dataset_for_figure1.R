# Clean environment
rm(list = ls())

# Libraries
library(dplyr)

# Set working directory
setwd("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data")

# Load dataset
dataset <- readRDS("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/data3_20250514.rds")

# Read in the list of IDs
id_list <- read.table("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/merged_ids.rel.id", header = FALSE)$V1

# Subset the dataset
dataset_subset <- dataset %>%
  filter(ID %in% id_list) %>%     # ensure 'ID' matches your dataset
  select(ID, DP0a, SP0a, PP0a, ethn1_consolidated)

# Save as CSV
write.csv(dataset_subset, "subset_dataset_figure1.csv", row.names = FALSE)
