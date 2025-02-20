rm(list=ls()); gc()
set.seed(1123)

# Load necessary libraries
library(dplyr)
library(readr)

# Set base directory
base_dir <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/"

# Load the scaled dataset
load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/scaled_dataset_20250106.Rdata")

# Initialize an empty data frame to store results
results <- data.frame(Ethnicity = character(), Trait = character(), Model = character(), R_squared = numeric(), Correlation = numeric(), stringsAsFactors = FALSE)

# Traits to process
traits <- c("SP", "DP", "PP")
# Ethnicities to process
ethnicities <- c("asian", "mixed", "black", "white", "chinese")
# Models to process
models <- c("X1", "X1_X2", "X1_X2_E", "X1_X2_G", "X1_X2_G_E", "X1_X2_G_E_GE", "X1_X2_G_E_GE_GEselect")

# Loop over ethnicities, traits, and models
for (ethn in ethnicities) {
  for (trait in traits) {
    for (model in models) {
      
      # Define file paths
      preds_file <- file.path(base_dir, paste0("PREDs_", trait, "_ethn_", ethn, "_", model, ".csv"))
      ethn_file <- file.path(base_dir, paste0(ethn, "_ids.txt"))
      
      # Load prediction dataset
      if (!file.exists(preds_file) || !file.exists(ethn_file)) {
        cat("Skipping Ethnicity", ethn, "Trait", trait, "Model", model, "due to missing files\n")
        next
      }
      
      preds_df <- read_csv(preds_file)
      
      # Load ethnicity IDs
      ethn_ids <- read.table(ethn_file)$V1
      
      # Subset the dataset based on Ethnicity IDs
      preds_subset <- preds_df %>%
        filter(ID %in% ethn_ids)
      
      # Fill missing values from the scaled dataset
      preds_filled <- preds_subset %>%
        left_join(dataset %>% select(ID, paste0(trait, "0s")), by = "ID") %>%
        mutate(Observed = ifelse(is.na(Observed), get(paste0(trait, "0s")), Observed)) %>%
        select(-matches(paste0(trait, "0s"))) # Remove extra columns after filling the missing values
      
      # Save the updated dataset
      write_csv(preds_filled, file.path(base_dir, paste0("Updated_PREDs_", trait, "_ethn_", ethn, "_", model, ".csv")))
      
      # Compute R^2 and correlation
      r_squared <- cor(preds_filled$Observed, preds_filled$Predicted, use = "complete.obs")^2
      correlation <- cor(preds_filled$Observed, preds_filled$Predicted, use = "complete.obs")
      
      # Append results
      results <- rbind(results, data.frame(Ethnicity = ethn, Trait = trait, Model = model, R_squared = r_squared, Correlation = correlation))
      
      # Print values for logging
      cat("Ethnicity:", ethn, "Trait:", trait, "Model:", model, "R^2:", r_squared, "Correlation:", correlation, "\n")
    }
  }
}

# Save results to CSV
write_csv(results, file.path(base_dir, "all_model_summary_ethnicity.csv"))

cat("Results saved to model_performance_summary_ethnicity.csv\n")
