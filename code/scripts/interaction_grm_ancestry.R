#Interaction GRM (One ancestry, All Environments) - Dynamic BP Handling ##

rm(list=ls()); gc()
set.seed(1123)

# Load necessary libraries
library(matrixcalc)
library(genio)
library(data.table)
library(Matrix)
library(optparse)

# Command-line options
option_list <- list(
  make_option(c("-d", "--dir"), type = "character", default = NULL,
              help = "path to the working directory", metavar = "character"),
  make_option(c("-b", "--bfile"), type = "character", default = NULL,
              help = "Path to the plink geno file", metavar = "character"),
  make_option(c("-e", "--envmat"), type = "character", default = NULL,
              help = "path to the env matrix", metavar = "character"),
  make_option(c("-s", "--snps"), type = "character", default = NULL,
              help = "path to filtered snps", metavar = "character"),
  make_option(c("-f", "--ancestry"), type = "character", default = NULL,
              help = "ancestry num", metavar = "character"),
  make_option(c("-p", "--bp"), type = "character", default = NULL,
              help = "BP", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output directory for saving results", metavar = "character")
  
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

print(opt)

# Define list of known interaction terms (environment terms)
interaction_terms <- c("Townsend", "act0_d", "TVtime", "sleep_d", "smoking_now", "veg_cook", "fish_oily",
                       "fish_lean", "meat_proc", "poultry", "beef", "lamb", "pork", "cheese", "salt", "tea",
                       "alc1", "waist", "getup", "coffee", "smoked_past", "BFP", "sleep_dev")

# Define parameters from command-line arguments
ancestry <- opt$ancestry
bp_type <- opt$bp
plink_results_dir <- opt$dir
geno_data_path <- opt$bfile
env_data_path <- opt$envmat
output_dir <- opt$output  # Ensure the output directory is correctly used

# Load PLINK data
plink_data <- read_plink(geno_data_path)
snp_matrix <- plink_data$X  # Extract genotype matrix
snp_names <- plink_data$bim$id  # SNP names

# Load environmental data
load(env_data_path)
env_matrix <- Emat

# Convert space-separated SNP files into a list
snp_files_passed <- unlist(strsplit(opt$snps, " "))

# Get dynamically found filtered files
filtered_files <- list.files(path = plink_results_dir, pattern = paste0("_", ancestry, ".", bp_type, "0s_filtered.txt$"), full.names = TRUE)

# Combine both lists (deduplicating if necessary)
filtered_files <- unique(c(filtered_files, snp_files_passed))  # ✅ Ensure all provided files are included
print(paste("📂 Using filtered SNP files:", paste(filtered_files, collapse=", ")))  # Debugging step

# Initialize final GRM for this ancestry and BP
num_individuals <- ncol(snp_matrix)
K_ancestry_grm <- matrix(0, num_individuals, num_individuals)  # Pre-allocate memory
num_envs_used <- 0  # Counter for environments used

# Loop through each known environment term
for (env_term in interaction_terms) {
  
  print(paste("🔹 Processing environment:", env_term, "for ancestry:", ancestry, "BP:", bp_type))
  
  # Identify the corresponding filtered file for this ancestry, BP, and environment
  filtered_file <- filtered_files[grepl(paste0("gxe_", env_term, "_", ancestry, ".", bp_type, "0s_filtered.txt$"), filtered_files)]
  
  if (length(filtered_file) == 0) {
    print(paste("⚠️ No filtered SNP file found for:", env_term, "in ancestry", ancestry, "BP:", bp_type))
    next
  }
  
  # Load SNPs from the filtered file
  gxe_results <- fread(filtered_file)
  env_snps <- gxe_results$ID  # Extract SNP column
  
  # Filter SNP matrix to keep only relevant SNPs
  snp_indices <- which(snp_names %in% env_snps)
  if (length(snp_indices) == 0) {
    print(paste("⚠️ No SNPs found in genotype data for:", env_term, "in ancestry", ancestry, "BP:", bp_type))
    next
  }
  
  snp_matrix_filtered <- snp_matrix[snp_indices, , drop=FALSE]
  
  # Extract the corresponding environmental variable
  if (!env_term %in% colnames(env_matrix)) {
    print(paste("⚠️ Skipping:", env_term, "- not found in environmental matrix"))
    next
  }
  
  env_vector <- env_matrix[, env_term]
  
  # Compute Environment GRM
  compute_env_grm_matrix <- function(env_vector) {
    outer(env_vector, env_vector) / length(env_vector)
  }
  
  K_env <- compute_env_grm_matrix(env_vector)
  
  # Pre-allocate memory for the SNP × ENV GRM for this environment
  K_snp_env <- matrix(0, num_individuals, num_individuals)
  
  # Compute SNP × ENV GRM iteratively
  compute_snp_grm_and_add <- function(snp_matrix, K_env, K_snp_env) {
    num_snps <- nrow(snp_matrix)
    
    for (s in 1:num_snps) {
      # Compute allele frequency
      p <- mean(snp_matrix[s, ], na.rm = TRUE) / 2
      
      # Center the SNP genotypes
      centered <- snp_matrix[s, ] - (2 * p)
      centered[is.na(centered)] <- 0
      
      # Compute SNP GRM for SNP s
      K_snp <- outer(centered, centered) / (2 * p * (1 - p))
      
      # Compute Hadamard product using `hadamard()`
      K_snp_env <- K_snp_env + hadamard.prod(K_snp, K_env)  # ✅ Using `hadamard()`
    }
    
    return(K_snp_env / num_snps)  # Averaging over SNPs
  }
  
  # Compute SNP × ENV GRM for this environment
  K_snp_env <- compute_snp_grm_and_add(snp_matrix_filtered, K_env, K_snp_env)
  
  # Accumulate into the full ancestry GRM
  K_ancestry_grm <- K_ancestry_grm + K_snp_env
  num_envs_used <- num_envs_used + 1
}

# Compute the final full averaged GRM for this ancestry and BP
if (num_envs_used > 0) {
  K_ancestry_grm <- K_ancestry_grm / num_envs_used  # ✅ Averaging over environments
  ancestry_output_file <- file.path(output_dir, paste0("Hadamard_GRM_ancestry", ancestry, "_", bp_type, ".RData"))
  save(K_ancestry_grm, file = ancestry_output_file)
  print(paste("💾 Saved Hadamard GRM for ancestry", ancestry, "BP:", bp_type, "to", ancestry_output_file))
} else {
  print(paste("⚠️ No environments were successfully processed for ancestry", ancestry, "BP:", bp_type, ". Full ancestry GRM was not computed."))
}
