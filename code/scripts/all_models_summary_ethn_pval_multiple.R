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
    # Load ethnicity IDs
    ethn_file <- file.path(base_dir, paste0(ancestry, "_ids.txt"))
    if (!file.exists(ethn_file)) {
      cat("Missing ancestry ID file for", ancestry, "\n")
      next
    }

    ethn_ids <- read.table(ethn_file)$V1  # Read ethnicity IDs correctly

    for (model in models) {
      if (model == "X1_X2_G_E_GE_GEselect") {
        for (pval in pvals) {
          preds_file <- file.path(base_dir,
                                  paste0("PREDs_", trait, "_ethn_", ancestry, "pval", pval, "_", model, ".csv"))

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

          preds_df <- read_csv(preds_file)

          # Keep only the individuals in ethn_ids
          preds_filtered <- preds_df %>% filter(ID %in% ethn_ids)

          # Fill missing values from the scaled dataset
          preds_filled <- preds_filtered %>%
            left_join(dataset %>% select(ID, paste0(trait, "0s")), by = "ID") %>%
            mutate(Observed = ifelse(is.na(Observed), get(paste0(trait, "0s")), Observed)) %>%
            select(-matches(paste0(trait, "0s")))

          # Save updated dataset
          updated_file <- file.path(base_dir,
                                    paste0("Updated_PREDs_", trait, "_ethn_", ancestry, "pval", pval, "_", model, ".csv"))
          write_csv(preds_filled, updated_file)

          # Compute R² and correlation
          r_squared <- cor(preds_filled$Observed, preds_filled$Predicted, use = "complete.obs")^2
          correlation <- cor(preds_filled$Observed, preds_filled$Predicted, use = "complete.obs")

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
        preds_file <- file.path(base_dir, paste0("PREDs_", trait, "_ethn_", ancestry, "_", model, ".csv"))

        if (!file.exists(preds_file)) {
          cat("Skipping Trait", trait, "Ancestry", ancestry, "Model", model, "due to missing file\n")
          results <- rbind(results, data.frame(Ancestry = ancestry,
                                               Trait = trait,
                                               Model = model,
                                               Pval = NA,
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
                                             Pval = NA,
                                             R_squared = r_squared,
                                             Correlation = correlation,
                                             stringsAsFactors = FALSE))

        cat("Trait:", trait, "Ancestry:", ancestry, "Model:", model,
            "R^2:", r_squared, "Correlation:", correlation, "\n")
      }
    }
  }
}

# Save the summary results
summary_file <- file.path(base_dir, "all_model_summary_pval_ethn_corrected.csv")
write_csv(results, summary_file)

cat("Results saved to", summary_file, "\n")
