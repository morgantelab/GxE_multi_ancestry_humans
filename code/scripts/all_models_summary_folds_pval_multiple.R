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
results <- data.frame(Fold = integer(),
                      Trait = character(),
                      Model = character(),
                      Pval = numeric(),
                      R_squared = numeric(),
                      Correlation = numeric(),
                      stringsAsFactors = FALSE)

# Define traits, folds, and models
traits <- c("SP", "DP", "PP")
folds <- 1:5
models <- c("X1", "X1_X2", "X1_X2_E", "X1_X2_G", "X1_X2_G_E", "X1_X2_G_E_GE", "X1_X2_G_E_GE_GEselect")
pvals <- c(1e-05, 1e-06, 0.001)  # Only applicable for GEselect

# Loop over folds, traits, and models
for (fold_num in folds) {
  for (trait in traits) {
    for (model in models) {
      
      if (model == "X1_X2_G_E_GE_GEselect") {
        # Loop over p-values only for GEselect
        for (pval in pvals) {
          # Define file paths with the pval component in the filename
          preds_file <- file.path(base_dir, paste0("PREDs_", trait, "_Fold_", fold_num, "pval", pval, "_", model, ".csv"))
          fold_file <- file.path(base_dir, paste0("Fold_", fold_num, ".rds"))
          
          # Check if required files exist
          if (!file.exists(preds_file) || !file.exists(fold_file)) {
            cat("Skipping Fold", fold_num, "Trait", trait, "Pval", pval, "Model", model, "due to missing files\n")
            results <- rbind(results, data.frame(Fold = fold_num,
                                                 Trait = trait,
                                                 Model = model,
                                                 Pval = pval,
                                                 R_squared = NA,
                                                 Correlation = NA,
                                                 stringsAsFactors = FALSE))
            next
          }
          
          # Load prediction dataset and fold IDs
          preds_df <- read_csv(preds_file)
          fold <- readRDS(fold_file)
          fold_ids <- fold$ID
          
          # Subset predictions based on fold IDs
          preds_subset <- preds_df %>% 
            filter(ID %in% fold_ids)
          
          # Fill missing values from the scaled dataset
          preds_filled <- preds_subset %>%
            left_join(dataset %>% select(ID, paste0(trait, "0s")), by = "ID") %>%
            mutate(Observed = ifelse(is.na(Observed), get(paste0(trait, "0s")), Observed)) %>%
            select(-matches(paste0(trait, "0s")))
          
          # Save the updated dataset with the pval included in the filename
          updated_file <- file.path(base_dir, paste0("Updated_PREDs_", trait, "_Fold_", fold_num, "pval", pval, "_", model, ".csv"))
          write_csv(preds_filled, updated_file)
          
          # Compute R² and correlation
          r_squared <- cor(preds_filled$Observed, preds_filled$Predicted, use = "complete.obs")^2
          correlation <- cor(preds_filled$Observed, preds_filled$Predicted, use = "complete.obs")
          
          # Append results
          results <- rbind(results, data.frame(Fold = fold_num,
                                               Trait = trait,
                                               Model = model,
                                               Pval = pval,
                                               R_squared = r_squared,
                                               Correlation = correlation,
                                               stringsAsFactors = FALSE))
          
          cat("Fold:", fold_num, "Trait:", trait, "Pval:", pval, "Model:", model,
              "R^2:", r_squared, "Correlation:", correlation, "\n")
        }
      } else {
        # For all other models (do not include a pval component in the filename)
        current_pval <- NA
        preds_file <- file.path(base_dir, paste0("PREDs_", trait, "_Fold_", fold_num, "_", model, ".csv"))
        fold_file <- file.path(base_dir, paste0("Fold_", fold_num, ".rds"))
        
        if (!file.exists(preds_file) || !file.exists(fold_file)) {
          cat("Skipping Fold", fold_num, "Trait", trait, "Model", model, "due to missing files\n")
          results <- rbind(results, data.frame(Fold = fold_num,
                                               Trait = trait,
                                               Model = model,
                                               Pval = current_pval,
                                               R_squared = NA,
                                               Correlation = NA,
                                               stringsAsFactors = FALSE))
          next
        }
        
        # Load prediction dataset and fold IDs
        preds_df <- read_csv(preds_file)
        fold <- readRDS(fold_file)
        fold_ids <- fold$ID
        
        # Subset predictions based on fold IDs
        preds_subset <- preds_df %>%
          filter(ID %in% fold_ids)
        
        # Fill missing values from the scaled dataset
        preds_filled <- preds_subset %>%
          left_join(dataset %>% select(ID, paste0(trait, "0s")), by = "ID") %>%
          mutate(Observed = ifelse(is.na(Observed), get(paste0(trait, "0s")), Observed)) %>%
          select(-matches(paste0(trait, "0s")))
        
        # Save the updated dataset without the pval in the filename
        updated_file <- file.path(base_dir, paste0("Updated_PREDs_", trait, "_Fold_", fold_num, "_", model, ".csv"))
        write_csv(preds_filled, updated_file)
        
        # Compute R² and correlation
        r_squared <- cor(preds_filled$Observed, preds_filled$Predicted, use = "complete.obs")^2
        correlation <- cor(preds_filled$Observed, preds_filled$Predicted, use = "complete.obs")
        
        # Append results
        results <- rbind(results, data.frame(Fold = fold_num,
                                             Trait = trait,
                                             Model = model,
                                             Pval = current_pval,
                                             R_squared = r_squared,
                                             Correlation = correlation,
                                             stringsAsFactors = FALSE))
        
        cat("Fold:", fold_num, "Trait:", trait, "Model:", model,
            "R^2:", r_squared, "Correlation:", correlation, "\n")
      }
    }
  }
}

# Save the summary results (missing values are represented as NA)
summary_file <- file.path(base_dir, "all_model_summary_folds_pval.csv")
write_csv(results, summary_file)

cat("Results saved to", summary_file, "\n")
