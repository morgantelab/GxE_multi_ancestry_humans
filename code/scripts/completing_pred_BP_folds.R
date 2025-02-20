rm(list=ls()); gc()
set.seed(1123)

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)

setwd("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model")

# Step 1: Load the datasets
# Load the prediction dataset
preds_df <- read_csv("PREDs_SP_Fold_3_X1_X2_G_E_GE_GEselect.csv")

# Load the IDs from Fold_1
fold <- readRDS("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/Fold_3.rds")
fold_ids <- fold$ID

# Load the scaled dataset
load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/scaled_dataset_20250106.Rdata")

# Step 2: Subset the dataset based on Fold_5 IDs
preds_subset <- preds_df %>%
  filter(ID %in% fold_ids)

# Step 3: Fill missing values from the scaled dataset
# Subset `preds_subset` IDs from the `dataset` and fill the `Observed` column
preds_filled <- preds_subset %>%
  left_join(dataset %>% select(ID, SP0s), by = "ID") %>%
  mutate(Observed = ifelse(is.na(Observed), SP0s, Observed)) %>%
  select(-SP0s) # Remove SP0s column after filling the missing values

# Step 4: Save the updated dataset
write_csv(preds_filled, "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/Updated_PREDs_SP_Fold_2_X1_X2_G_E_GE_GEselect.csv")

# Step 5: Calculate R^2 and correlation
# Assuming the dataset has columns `Observed` and `Predicted`
r_squared <- cor(preds_filled$Observed, preds_filled$Predicted, use = "complete.obs")^2
correlation <- cor(preds_filled$Observed, preds_filled$Predicted, use = "complete.obs")

# Print the values
cat("R^2:", r_squared, "\n")
cat("Correlation:", correlation, "\n")
