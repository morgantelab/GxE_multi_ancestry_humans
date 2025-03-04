# Clear workspace and set seed for reproducibility
rm(list=ls()); gc()
set.seed(1123)

# Load required packages
library(dplyr)
library(purrr)
library(stringr)

# Define parameters
traits <- c("SP", "DP", "PP")
ancestries <- c("asian", "mixed", "black", "chinese", "white")
models <- c("X1", "X1_X2", "X1_X2_E", "X1_X2_G", "X1_X2_G_E", "X1_X2_G_E_GE", "X1_X2_G_E_GE_GEselect")
pvals <- c("1e-05", "1e-06", "1e-03")  # Ensure pvals are character strings
folds <- 1:5  # Folds 1 to 5

# Define directories
save_dir <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model"
data_dir <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model"

# Load the scaled dataset
load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/scaled_dataset_20250106.Rdata")

# Ensure dataset is loaded correctly
if (!exists("dataset")) {
  stop("Error: dataset not found after loading the Rdata file.")
}

# Function to generate filenames
generate_filenames <- function(trait, ancestry, model, fold = NULL) {
  if (!is.null(fold)) {  # For fold-specific files
    if (model == "X1_X2_G_E_GE_GEselect") {
      return(sprintf("PREDs_%s_Fold_%d_pval_%s_%s.csv", trait, fold, pvals, model))
    } else {
      return(sprintf("PREDs_%s_Fold_%d_%s.csv", trait, fold, model))
    }
  } else {  # For non-fold files
    if (model == "X1_X2_G_E_GE_GEselect") {
      return(sprintf("PREDs_%s_ethn_%s_pval_%s_%s.csv", trait, ancestry, pvals, model))
    } else {
      return(sprintf("PREDs_%s_ethn_%s_%s.csv", trait, ancestry, model))
    }
  }
}

# Generate filenames for ancestry-based files
file_combinations_ancestry <- expand.grid(traits, ancestries, models, stringsAsFactors = FALSE) %>%
  pmap(~ generate_filenames(..1, ..2, ..3)) %>%
  unlist()

# Generate filenames for fold-based files
file_combinations_folds <- expand.grid(traits, folds, models, stringsAsFactors = FALSE) %>%
  pmap(~ generate_filenames(..1, NULL, ..3, ..2)) %>%
  unlist()

# Combine all filenames
file_names <- c(file_combinations_ancestry, file_combinations_folds)

# Construct full file paths
file_paths <- file.path(data_dir, file_names)

# Function to read a file and return both data and summary
read_and_summarize <- function(file_path) {
  if (file.exists(file_path)) {
    data <- read.csv(file_path)
    message(sprintf("Loaded: %s (%d rows, %d columns)", basename(file_path), nrow(data), ncol(data)))
    return(list(file = basename(file_path), data = data, rows = nrow(data), columns = ncol(data)))
  } else {
    message(sprintf("File not found: %s", basename(file_path)))
    return(NULL)
  }
}

# Read files and store results
results <- lapply(file_paths, read_and_summarize)

# Filter out NULLs
loaded_files <- results[!sapply(results, is.null)]

# Extract data and summaries separately
preds_list <- setNames(lapply(loaded_files, function(x) x$data), sapply(loaded_files, function(x) x$file))
summary_list <- do.call(rbind, lapply(loaded_files, function(x) data.frame(File = x$file, Rows = x$rows, Columns = x$columns)))

# Display summary of loaded files
print(summary_list)

# Function to filter out rows where 'Observed' is NA
filter_observed_na <- function(df) {
  if ("Observed" %in% colnames(df)) {
    return(df %>% filter(is.na(Observed)))  # Keep only rows where Observed is NA
  } else {
    warning("Observed column not found in dataset.")
    return(df)  # Return unmodified if column is missing
  }
}

# Apply the filtering to all files in preds_list
filtered_preds_list <- lapply(preds_list, filter_observed_na)

# Function to extract trait from filename
extract_trait <- function(filename) {
  trait_match <- str_match(filename, "PREDs_([A-Za-z]+)_")[,2]
  return(trait_match)
}

# Function to calculate R-squared
calculate_r_squared <- function(observed, predicted) {
  valid_indices <- !is.na(observed) & !is.na(predicted)
  if (sum(valid_indices) > 1) {
    return(cor(observed[valid_indices], predicted[valid_indices])^2)
  } else {
    return(NA)  # Return NA if not enough valid data points
  }
}

# Function to replace NA in Observed column using dataset and compute R-squared
update_observed_values <- function(df, filename) {
  if (!"ID" %in% colnames(df) || !"Observed" %in% colnames(df) || !"Predicted" %in% colnames(df)) {
    warning(sprintf("Skipping %s: Missing ID, Observed, or Predicted column.", filename))
    return(list(data = df, r_squared = NA))
  }
  
  # Extract trait from filename
  trait <- extract_trait(filename)
  trait_column <- paste0(trait, "0s")  # Match dataset column (e.g., DP → DP0s)
  
  # Check if the trait column exists in dataset
  if (!(trait_column %in% colnames(dataset))) {
    warning(sprintf("Skipping %s: Trait column %s not found in dataset.", filename, trait_column))
    return(list(data = df, r_squared = NA))
  }
  
  # Replace NA values in Observed column where ID matches
  id_match <- match(df$ID, dataset$ID)  # Find matching indices
  valid_indices <- which(!is.na(id_match) & is.na(df$Observed))  # Only replace NA values
  
  if (length(valid_indices) > 0) {
    df$Observed[valid_indices] <- dataset[[trait_column]][id_match[valid_indices]]
  }
  
  # Compute R-squared after updating Observed
  r_squared <- calculate_r_squared(df$Observed, df$Predicted)
  
  return(list(data = df, r_squared = r_squared))
}

# Apply NA replacement using dataset and compute R-squared
updated_results <- mapply(update_observed_values, filtered_preds_list, names(filtered_preds_list), SIMPLIFY = FALSE)

# Extract updated datasets and R-squared values
updated_preds_list <- lapply(updated_results, function(x) x$data)
r_squared_values <- sapply(updated_results, function(x) x$r_squared)

# Save updated datasets
for (file_name in names(updated_preds_list)) {
  new_file_name <- paste0("Updated_predictions_", file_name)
  new_file_path <- file.path(save_dir, new_file_name)
  write.csv(updated_preds_list[[file_name]], new_file_path, row.names = FALSE)
  message(sprintf("Saved updated file: %s", new_file_name))
}

# Create and save final summary file including R-squared
updated_summary <- data.frame(
  File = names(updated_preds_list),
  Remaining_Rows = sapply(updated_preds_list, nrow),
  R_Squared = r_squared_values
)

summary_file_path <- file.path(save_dir, "Updated_predictions_summary.csv")
write.csv(updated_summary, summary_file_path, row.names = FALSE)
message("Saved summary file: Updated_predictions_summary.csv")

# Print final summary
print(updated_summary)
