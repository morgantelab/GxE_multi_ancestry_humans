rm(list=ls()); gc()
set.seed(1123)

# Load necessary libraries
library(dplyr)
library(readr)

# Set base directory
base_dir <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model"

# Load the scaled dataset
load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/scaled_dataset_20250106.Rdata")

# Initialize an empty data frame to store results
results <- data.frame(Fold = integer(), Trait = character(), R_squared = numeric(), Correlation = numeric(), stringsAsFactors = FALSE)

# Traits to process
traits <- c("SP", "DP", "PP")

# Loop over folds 1 to 5
for (fold_num in 1:5) {
  for (trait in traits) {
    
    # Define file paths
    preds_file <- file.path(base_dir, paste0("PREDs_", trait, "_Fold_", fold_num, "_X1_X2_G_E_GE_GEselect.csv"))
    fold_file <- file.path(base_dir, paste0("Fold_", fold_num, ".rds"))
    
    # Load prediction dataset
    if (!file.exists(preds_file) || !file.exists(fold_file)) {
      cat("Skipping Fold", fold_num, "for trait", trait, "due to missing files\n")
      next
    }
    
    preds_df <- read_csv(preds_file)
    
    # Load fold IDs
    fold <- readRDS(fold_file)
    fold_ids <- fold$ID
    
    # Subset the dataset based on Fold IDs
    preds_subset <- preds_df %>%
      filter(ID %in% fold_ids)
    
    # Fill missing values from the scaled dataset
    preds_filled <- preds_subset %>%
      left_join(dataset %>% select(ID, paste0(trait, "0s")), by = "ID") %>%
      mutate(Observed = ifelse(is.na(Observed), get(paste0(trait, "0s")), Observed)) %>%
      select(-matches(paste0(trait, "0s"))) # Remove extra columns after filling the missing values
    
    # Save the updated dataset
    write_csv(preds_filled, file.path(base_dir, paste0("Updated_PREDs_", trait, "_Fold_", fold_num, "_X1_X2_G_E_GE_GEselect.csv")))
    
    # Compute R^2 and correlation
    r_squared <- cor(preds_filled$Observed, preds_filled$Predicted, use = "complete.obs")^2
    correlation <- cor(preds_filled$Observed, preds_filled$Predicted, use = "complete.obs")
    
    # Append results
    results <- rbind(results, data.frame(Fold = fold_num, Trait = trait, R_squared = r_squared, Correlation = correlation))
    
    # Print values for logging
    cat("Fold:", fold_num, "Trait:", trait, "R^2:", r_squared, "Correlation:", correlation, "\n")
  }
}

# Save results to CSV
write_csv(results, file.path(base_dir, "model_performance_summary_folds.csv"))

cat("Results saved to model_performance_summary.csv\n")
