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
results <- data.frame(Ancestry = character(),
                      Trait = character(),
                      Model = character(),
                      Pval = numeric(),
                      R_squared = numeric(),
                      Correlation = numeric(),
                      stringsAsFactors = FALSE)

# Define traits, ancestries, models, and p-value thresholds
traits <- c("SP", "DP", "PP")
ancestries <- c("asian", "mixed", "black", "chinese", "white")
models <- c("X1", "X1_X2", "X1_X2_E", "X1_X2_G", "X1_X2_G_E", "X1_X2_G_E_GE", "X1_X2_G_E_GE_GEselect")
pvals <- c(1e-05, 1e-06, 0.001)  # Only applicable for GEselect

# Loop over traits, ancestries, and models
for (trait in traits) {
  for (ancestry in ancestries) {
    for (model in models) {
      
      if (model == "X1_X2_G_E_GE_GEselect") {
        # Loop over p-values only for GEselect model
        for (pval in pvals) {
          # Construct file name including the pval component
          preds_file <- file.path(base_dir, 
                                  paste0("PREDs_", trait, "_ethn_", ancestry, "pval", pval, "_", model, ".csv"))
          
          # Check if the prediction file exists; if not, record NA values in results
          if (!file.exists(preds_file)) {
            cat("Skipping Trait", trait, "Ancestry", ancestry, "Pval", pval, "Model", model, "due to missing file\n")
            results <- rbind(results, data.frame(Ancestry = ancestry,
                                                 Trait = trait,
                                                 Model = model,
                                                 Pval = pval,
                                                 R_squared = NA,
                                                 Correlation = NA,
                                                 stringsAsFactors = FALSE))
            next
          }
          
          # Load prediction dataset
          preds_df <- read_csv(preds_file)
          
          # Fill missing values from the scaled dataset
          preds_filled <- preds_df %>%
            left_join(dataset %>% select(ID, paste0(trait, "0s")), by = "ID") %>%
            mutate(Observed = ifelse(is.na(Observed), get(paste0(trait, "0s")), Observed)) %>%
            select(-matches(paste0(trait, "0s")))
          
          # Save the updated dataset with the pval in the filename
          updated_file <- file.path(base_dir, 
                                    paste0("Updated_PREDs_", trait, "_ethn_", ancestry, "pval", pval, "_", model, ".csv"))
          write_csv(preds_filled, updated_file)
          
          # Compute R² and correlation
          r_squared <- cor(preds_filled$Observed, preds_filled$Predicted, use = "complete.obs")^2
          correlation <- cor(preds_filled$Observed, preds_filled$Predicted, use = "complete.obs")
          
          # Append results
          results <- rbind(results, data.frame(Ancestry = ancestry,
                                               Trait = trait,
                                               Model = model,
                                               Pval = pval,
                                               R_squared = r_squared,
                                               Correlation = correlation,
                                               stringsAsFactors = FALSE))
          
          cat("Trait:", trait, "Ancestry:", ancestry, "Pval:", pval, "Model:", model,
              "R^2:", r_squared, "Correlation:", correlation, "\n")
        }
      } else {
        # For non-GEselect models, no pval component is in the file name
        current_pval <- NA
        preds_file <- file.path(base_dir, paste0("PREDs_", trait, "_ethn_", ancestry, "_", model, ".csv"))
        
        if (!file.exists(preds_file)) {
          cat("Skipping Trait", trait, "Ancestry", ancestry, "Model", model, "due to missing file\n")
          results <- rbind(results, data.frame(Ancestry = ancestry,
                                               Trait = trait,
                                               Model = model,
                                               Pval = current_pval,
                                               R_squared = NA,
                                               Correlation = NA,
                                               stringsAsFactors = FALSE))
          next
        }
        
        # Load prediction dataset
        preds_df <- read_csv(preds_file)
        
        # Fill missing values from the scaled dataset
        preds_filled <- preds_df %>%
          left_join(dataset %>% select(ID, paste0(trait, "0s")), by = "ID") %>%
          mutate(Observed = ifelse(is.na(Observed), get(paste0(trait, "0s")), Observed)) %>%
          select(-matches(paste0(trait, "0s")))
        
        # Save the updated dataset without the pval component in the filename
        updated_file <- file.path(base_dir, paste0("Updated_PREDs_", trait, "_ethn_", ancestry, "_", model, ".csv"))
        write_csv(preds_filled, updated_file)
        
        # Compute R² and correlation
        r_squared <- cor(preds_filled$Observed, preds_filled$Predicted, use = "complete.obs")^2
        correlation <- cor(preds_filled$Observed, preds_filled$Predicted, use = "complete.obs")
        
        # Append results
        results <- rbind(results, data.frame(Ancestry = ancestry,
                                             Trait = trait,
                                             Model = model,
                                             Pval = current_pval,
                                             R_squared = r_squared,
                                             Correlation = correlation,
                                             stringsAsFactors = FALSE))
        
        cat("Trait:", trait, "Ancestry:", ancestry, "Model:", model,
            "R^2:", r_squared, "Correlation:", correlation, "\n")
      }
    }
  }
}

# Save the summary results (missing values are represented as NA)
summary_file <- file.path(base_dir, "all_model_summary_folds_pval_ethn.csv")
write_csv(results, summary_file)

cat("Results saved to", summary_file, "\n")
