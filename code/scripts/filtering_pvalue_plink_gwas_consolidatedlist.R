rm(list=ls()); gc()
set.seed(1123)

# Load necessary library
library(data.table)

# Define directory where filtered files are stored
filtered_dir <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/gwas_snp_env/"

# List all filtered result files
filtered_files <- list.files(path = filtered_dir, pattern = "\\.filtered\\.linear$", full.names = TRUE)

# Function to extract details from filename and count rows
extract_info <- function(file) {
  df <- fread(file, header = TRUE, select = 1)  # Read only one column for speed
  
  # Extract details from filename
  file_basename <- basename(file)
  parts <- strsplit(file_basename, "_")[[1]]
  
  # Identify parts of the filename (adjust based on your naming convention)
  trait <- parts[2]  # Trait name
  env <- parts[3]  # Environmental variable
  fold_ancestry <- gsub("\\.filtered\\.linear$", "", paste(parts[4:length(parts)], collapse = "_"))  # Fold or ancestry
  
  # Return structured data
  data.frame(Trait = trait, Env = env, `Fold/ancestry` = fold_ancestry, nrow = nrow(df))
}

# Apply function to all filtered files
summary_df <- do.call(rbind, lapply(filtered_files, extract_info))

# Sort summary by nrow (biggest to smallest)
summary_df <- summary_df[order(-summary_df$nrow), ]

# Display the summary
print(summary_df)

# Save summary as CSV
summary_file <- file.path(filtered_dir, "filtered_summary.csv")
fwrite(summary_df, summary_file, sep = ",")

cat("Summary saved to:", summary_file, "\n")
