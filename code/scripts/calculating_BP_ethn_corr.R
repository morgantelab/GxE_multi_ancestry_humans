# Clear the workspace
rm(list=ls()); gc()
# Load necessary libraries
library(dplyr)
library(readr)
library(tidyr)
# Set working directory
setwd("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model")
# Define parameters
models <- c("X1", "X1_X2", "X1_X2_E", "X1_X2_G", "X1_X2_G_E", "X1_X2_G_E_GE")
bps <- c("DP", "SP", "PP")
ethnicities <- c("White", "Black", "Asian", "Mixed", "Chinese")
# Create an empty dataframe for results
results <- data.frame(Model = character(),
                      BP = character(),
                      Ethnicity = character(),
                      R2 = numeric(),
                      Corr = numeric(),
                      stringsAsFactors = FALSE)
# Function to calculate R^2 and correlation
calculate_metrics <- function(model, bp, ethnicity) {
  # Load the prediction dataset
  preds_file <- sprintf("PREDs_%s_ethn_%s_%s.csv", bp, tolower(ethnicity), model)
  preds_df <- read_csv(preds_file)
  # Load the IDs
  ids_file <- sprintf("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/%s_ids.txt", tolower(ethnicity))
  ethn_ids <- read.table(ids_file)$V1
  # Subset the dataset
  preds_subset <- preds_df %>%
    filter(ID %in% ethn_ids)
  # Load the scaled dataset
  load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/scaled_dataset_20241025.Rdata")
  # Fill missing Observed values
  preds_filled <- preds_subset %>%
    left_join(dataset %>% select(ID, DP0s), by = "ID") %>%
    mutate(Observed = ifelse(is.na(Observed), DP0s, Observed)) %>%
    select(-DP0s)
  # Calculate metrics
  r_squared <- cor(preds_filled$Observed, preds_filled$Predicted, use = "complete.obs")^2
  correlation <- cor(preds_filled$Observed, preds_filled$Predicted, use = "complete.obs")
  return(list(R2 = r_squared, Corr = correlation))
}
# Loop through each combination of model, BP, and ethnicity
for (model in models) {
  for (bp in bps) {
    for (ethnicity in ethnicities) {
      metrics <- calculate_metrics(model, bp, ethnicity)
      results <- results %>%
        add_row(Model = model, BP = bp, Ethnicity = ethnicity,
                R2 = metrics$R2, Corr = metrics$Corr)
    }
  }
}
# Save the results to a CSV file
write_csv(results, "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/predictions_calcs/R2_Corr_Results_ethnicities.csv")