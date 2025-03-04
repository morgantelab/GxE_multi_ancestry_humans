rm(list=ls()); gc()
set.seed(1123)

# Load necessary libraries
library(dplyr)
library(readr)

# Load the scaled dataset
load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/scaled_dataset_20250106.Rdata")

# Define function to check ethnicity classification
check_ethnicity <- function(file_path, ethnicity_label, dataset) {
  # Read IDs from the first column of the text file and ensure they are character type
  ids <- read_delim(file_path, delim = " ", col_names = FALSE, col_types = cols(X1 = col_character()))$X1
  
  # Ensure dataset ID column is also character for accurate matching
  dataset <- dataset %>% mutate(ID = as.character(ID))
  
  # Filter dataset to keep only those IDs
  filtered_data <- dataset %>% filter(ID %in% ids)
  
  # Count how many match the expected ethnicity
  match_count <- sum(filtered_data$ethn1_consolidated == ethnicity_label, na.rm = TRUE)
  total_count <- length(ids)
  
  # Print summary
  cat(sprintf("Ethnicity Check: %s\n", ethnicity_label))
  cat(sprintf("Total IDs in %s: %d\n", file_path, total_count))
  cat(sprintf("Matched IDs: %d (%.2f%%)\n", match_count, (match_count / total_count) * 100))
  cat("\n")
  
  # Return a summary as a dataframe
  return(data.frame(
    Ethnicity = ethnicity_label,
    Total_IDs = total_count,
    Matched_IDs = match_count,
    Match_Percentage = ifelse(total_count > 0, (match_count / total_count) * 100, NA)
  ))
}

# Base directory where ID files are stored
base_dir <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/"

# List of ethnicity files and corresponding labels
ethnicity_files <- list(
  "white_ids.txt" = "White",
  "black_ids.txt" = "Black",
  "asian_ids.txt" = "Asian",
  "mixed_ids.txt" = "Mixed",
  "chinese_ids.txt" = "Chinese"
)

# Run checks for each ethnicity file
results <- lapply(names(ethnicity_files), function(file) {
  file_path <- file.path(base_dir, file)
  check_ethnicity(file_path, ethnicity_files[[file]], dataset)
})

# Combine results into a summary dataframe
summary_df <- bind_rows(results)

# Print the summary dataframe
print(summary_df)
