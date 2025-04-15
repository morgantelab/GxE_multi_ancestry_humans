rm(list=ls()); gc()
set.seed(1123)
library(data.table)

# Directory containing interaction term files
results_dir <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/gwas_snp_env/"

# List all interaction term CSV files
interaction_files <- list.files(results_dir, pattern="interaction_terms.*\\.csv$", full.names=TRUE)

# Loop through each interaction term file
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

cat("All files have been processed and filtered for lowest 10 P-values.\n")
