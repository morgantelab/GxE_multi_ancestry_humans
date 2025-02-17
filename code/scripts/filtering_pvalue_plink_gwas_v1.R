rm(list=ls()); gc()
set.seed(1123)

# Load necessary library
library(data.table)

# Define input and output directories
input_dir <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/gwas_snp_env/"
output_dir <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/gwas_snp_env/"

# List all GWAS result files (removes double slashes)
all_files <- list.files(path = input_dir, pattern = "^gxe_.*\\.glm\\.linear$", full.names = TRUE)
all_files <- gsub("//", "/", all_files)  # Remove any unwanted double slashes

# Function to process and filter a single GWAS file
filter_gwas <- function(file) {
  cat("Processing:", file, "\n")
  
  # Read the GWAS results
  df <- fread(file, header = TRUE)
  
  # Identify the P-value column dynamically
  p_col <- grep("P$", colnames(df), value = TRUE)
  if (length(p_col) == 0) {
    cat("Warning: No P-value column found in", file, "\n")
    return(NULL)
  }
  
  # Extract the environmental variable from the filename
  file_basename <- basename(file)
  env_var <- sub("gxe_([^_]+)_.*", "\\1", file_basename)  # Extract env name after 'gxe_'
  
  # Filter for interaction terms ADDx{env} and P-value < 1e-5
  df_filtered <- df[grepl(paste0("ADDx", env_var, "$"), df$TEST) & df[[p_col]] < 1e-5, ]
  
  if (nrow(df_filtered) == 0) {
    cat("No significant interactions found for", env_var, "in", file, "\n")
    return(NULL)
  }
  
  # Define correct output filename
  output_file <- file.path(output_dir, paste0(file_basename, ".filtered.linear"))
  
  # Save filtered results
  fwrite(df_filtered, output_file, sep = "\t")
  
  cat("Filtered results saved to:", output_file, "\n")
}

# Apply the filtering function to all files
lapply(all_files, filter_gwas)

cat("Batch filtering completed for all environments!\n")


