
rm(list=ls()); gc()
set.seed(1123)

# Load necessary libraries
library(dplyr)
library(readr)

# Load the scaled dataset
load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/scaled_dataset_20250106.Rdata")

# Load necessary libraries
library(dplyr)
library(readr)

# Define function to check ethnicity breakdown
ethnicity_breakdown <- function(file_path, dataset) {
  # Read only the first column as character
  ids <- read_delim(file_path, delim = "\t", col_names = FALSE, col_types = cols(X1 = col_character()))$X1
  
  # Ensure dataset ID column is also character for accurate matching
  dataset <- dataset %>% mutate(ID = as.character(ID))
  
  # Find matching IDs in dataset
  matched_data <- dataset %>% filter(ID %in% ids)
  
  # Check if any matches exist
  if (nrow(matched_data) == 0) {
    cat(sprintf("No matching IDs found in dataset for file: %s\n", file_path))
    return(data.frame(ethn1_consolidated = character(0), Count = integer(0)))
  }
  
  # Summarize counts of ethn1_consolidated
  ethnicity_summary <- matched_data %>%
    group_by(ethn1_consolidated) %>%
    summarise(Count = n(), .groups = "drop") %>%
    arrange(desc(Count))
  
  # Print summary
  cat(sprintf("Ethnicity breakdown for IDs in %s:\n", file_path))
  print(ethnicity_summary)
  
  # Return the summary as a dataframe
  return(ethnicity_summary)
}

# Path to Train_exclude_white.txt
file_path <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/gwas_snp_env/Train_exclude_white.txt"

# Run ethnicity breakdown check
ethnicity_summary_df <- ethnicity_breakdown(file_path, dataset)

# Display the summary dataframe
print(ethnicity_summary_df)
