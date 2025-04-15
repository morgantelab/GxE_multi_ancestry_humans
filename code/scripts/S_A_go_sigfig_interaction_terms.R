rm(list=ls()); gc()
set.seed(1123)
library(data.table)

# Directory containing interaction term files
results_dir <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/gwas_snp_env/"

# Environments of interest
envs <- c("sex", "age", "age2")

# Get all interaction term files
all_files <- list.files(results_dir, pattern="S_A_gxe_.*_interaction_terms\\.csv$", full.names=TRUE)

# Filter files based on <env> extracted from the filename
interaction_files <- all_files[
  grepl(paste0("S_A_gxe_(", paste(envs, collapse="|"), ")_full_dataset"), basename(all_files))
]

# Loop through each selected interaction term file
for (file in interaction_files) {
  cat("Processing file:", file, "\n")
  
  # Read file
  data <- fread(file)
  
  # Ensure P column is numeric
  data[, P := as.numeric(P)]
  
  # Remove rows with NA in P
  data <- data[!is.na(P)]
  
  # Sort by P-values and keep lowest 10
  data_sorted <- data[order(P)][1:min(.N, 10)]
  
  # Generate output filename
  output_filename <- gsub("interaction_terms", "interaction_terms_top10", basename(file))
  output_filepath <- file.path(results_dir, output_filename)
  
  # Save the filtered top 10 interactions
  fwrite(data_sorted, output_filepath)
  cat("Saved top 10 interaction terms to:", output_filepath, "\n")
}

cat("Filtered files for sex, age, and age2 have been processed.\n")
