rm(list=ls()); gc()
set.seed(1123)

# Load necessary library
library(data.table)

# Define input and output directories
input_dir <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/gwas_snp_env/"
output_dir <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/gwas_snp_env/"

# List all GWAS result files
all_files <- list.files(path = input_dir, pattern = "^gxe_.*\\.glm\\.linear$", full.names = TRUE)
all_files <- gsub("//", "/", all_files)  # Remove any unwanted double slashes

# Define the list of environmental variables
env_vars <- c("act0_d", "sleep_d", "smoking_now", "veg_cook", "fish_oily", "fish_lean", "meat_proc", "smoked_past", "sleep_dev")

# Function to filter GWAS results for a given environmental variable
filter_gwas <- function(env) {
  # Identify files relevant to this environmental variable
  env_files <- all_files[grepl(paste0("gxe_", env, "_"), all_files)]
  
  if (length(env_files) == 0) {
    cat("No GWAS files found for:", env, "\n")
    return(NULL)
  }
  
  for (file in env_files) {
    cat("Processing:", file, "\n")
    
    # Read GWAS results
    df <- fread(file, header = TRUE)
    
    # Identify the P-value column dynamically
    p_col <- grep("^P$", colnames(df), value = TRUE)
    if (length(p_col) == 0) {
      cat("Warning: No P-value column found in", file, "\n")
      next
    }
    
    # Filter for the specified environment term directly
    df_filtered <- df[TEST == paste0("ADDx", env) & df[[p_col]] < 1e-5, ]
    
    if (nrow(df_filtered) == 0) {
      cat("No significant interactions found for", env, "in", file, "\n")
      next
    }
    
    # Define output filename
    output_file <- file.path(output_dir, paste0(basename(file), "2.filtered.linear"))
    
    # Save filtered results
    fwrite(df_filtered, output_file, sep = "\t")
    
    cat("Filtered results saved to:", output_file, "\n")
  }
}

# Apply the filtering function to all environmental variables
lapply(env_vars, filter_gwas)

cat("Batch filtering completed for all environments!\n")
