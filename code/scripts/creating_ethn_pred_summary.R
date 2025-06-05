rm(list=ls()); gc()
set.seed(1123)

# Load necessary libraries
library(dplyr)
library(readr)

# Set base directory
base_dir <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/"
ancestry_dir <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/"

# Load the scaled dataset
dataset <- readRDS("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/data3_20250514.rds")

# Initialize an empty data frame to store results
results <- data.frame(Ancestry = character(),
                      Trait = character(),
                      Model = character(),
                      R_squared = numeric(),
                      Correlation = numeric(),
                      stringsAsFactors = FALSE)

# Define traits, ancestries, and models
traits <- c("SP", "DP", "PP")
ancestries <- c("asian", "mixed", "black", "chinese", "white")
models <- c("X1", "X1_X2", "X1_X2_E", "X1_X2_G", "X1_X2_G_E", "X1_X2_G_E_GE")

# Loop over traits, ancestries, and models
for (trait in traits) {
  for (ancestry in ancestries) {
    # Load ethnicity IDs
    ethn_file <- file.path(ancestry_dir, paste0(ancestry, "_ids.txt"))
    if (!file.exists(ethn_file)) {
      cat("Missing ancestry ID file for", ancestry, "\n")
      next
    }
    
    ethn_ids <- read.table(ethn_file)$V1  # Read ethnicity IDs correctly
    
    for (model in models) {
      preds_file <- file.path(base_dir, paste0("PREDs_", trait, "_ethn_", ancestry, "_", model, ".csv"))
      
      if (!file.exists(preds_file)) {
        cat("Skipping Trait", trait, "Ancestry", ancestry, "Model", model, "due to missing file\n")
        results <- rbind(results, data.frame(Ancestry = ancestry,
                                             Trait = trait,
                                             Model = model,
                                             R_squared = NA,
                                             Correlation = NA,
                                             stringsAsFactors = FALSE))
        next
      }
      
      preds_df <- read_csv(preds_file)
      
      # Keep only the individuals in ethn_ids
      preds_filtered <- preds_df %>% filter(ID %in% ethn_ids)
      
      # Fill missing values from the scaled dataset
      preds_filled <- preds_filtered %>%
        left_join(dataset %>% select(ID, paste0(trait, "0s")), by = "ID") %>%
        mutate(Observed = ifelse(is.na(Observed), get(paste0(trait, "0s")), Observed)) %>%
        select(-matches(paste0(trait, "0s")))
      
      # Save updated dataset
      updated_file <- file.path(base_dir, paste0("Updated_PREDs_", trait, "_ethn_", ancestry, "_", model, ".csv"))
      write_csv(preds_filled, updated_file)
      
      # Compute R² and correlation
      r_squared <- cor(preds_filled$Observed, preds_filled$Predicted, use = "complete.obs")^2
      correlation <- cor(preds_filled$Observed, preds_filled$Predicted, use = "complete.obs")
      
      results <- rbind(results, data.frame(Ancestry = ancestry,
                                           Trait = trait,
                                           Model = model,
                                           R_squared = r_squared,
                                           Correlation = correlation,
                                           stringsAsFactors = FALSE))
      
      cat("Trait:", trait, "Ancestry:", ancestry, "Model:", model,
          "R^2:", r_squared, "Correlation:", correlation, "\n")
    }
  }
}

# Save the summary results
summary_file <- file.path(base_dir, "all_model_summary_ethn_corrected.csv")
write_csv(results, summary_file)

cat("Results saved to", summary_file, "\n")
