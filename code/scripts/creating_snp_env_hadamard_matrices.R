# V5 Interaction GRM (Filtered SNPs for Waist ENV) - Memory Optimized ##

rm(list=ls()); gc()
set.seed(1123)

# Load required libraries
library(matrixcalc)
library(genio)
library(data.table)
library(Matrix)

# Load PLINK data using genio
plink_data <- read_plink("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/merged_geno_common_snps_selected_indiv")
snp_matrix <- plink_data$X  # Extract genotype matrix
snp_names <- plink_data$bim$id  # SNP names

# Load environment data
load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/Emat_20250106.RData")
env_matrix <- Emat

# Load list of interacting SNPs for waist environment
gxe_results <- fread("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/gwas_snp_env/gxe_waist_4.DP0s_filtered.txt")
waist_snps <- gxe_results$ID  # Extract SNP column

# Filter SNP matrix to keep only relevant SNPs
snp_indices <- which(snp_names %in% waist_snps)
if (length(snp_indices) == 0) {
  stop("Error: No SNPs found in snp_matrix that match waist_snps.")
}
snp_matrix_filtered <- snp_matrix[snp_indices, , drop=FALSE]  # Ensure matrix structure

# Select only the "waist" environment from the environment matrix
if (!"waist" %in% colnames(env_matrix)) {
  stop("Error: 'waist' environment variable not found in Emat.")
}
waist_env_vector <- env_matrix[, "waist"]  # Extract as a vector

# Function to compute Environment GRM
compute_env_grm_matrix <- function(env_vector) {
  outer(env_vector, env_vector) / length(env_vector)
}

# Compute Environment GRM for waist
K_env <- compute_env_grm_matrix(waist_env_vector)

# Pre-allocate memory for the final SNP × ENV GRM
num_individuals <- ncol(snp_matrix_filtered)
K_snp_env_final <- matrix(0, num_individuals, num_individuals)

# Compute SNP × ENV GRM iteratively, avoiding storage of intermediate matrices
compute_snp_grm_and_add <- function(snp_matrix, K_env, K_snp_env_final) {
  num_snps <- nrow(snp_matrix)

  for (s in 1:num_snps) {
    # Compute allele frequency
    p <- mean(snp_matrix[s, ], na.rm = TRUE) / 2

    # Center the SNP genotypes
    centered <- snp_matrix[s, ] - (2 * p)
    centered[is.na(centered)] <- 0  # Impute missing values

    # Compute SNP GRM on-the-fly
    K_snp <- outer(centered, centered) / (2 * p * (1 - p))

    # Add Hadamard product to final matrix
    K_snp_env_final <- K_snp_env_final + (K_snp * K_env)

  }

  # Average over SNPs
  return(K_snp_env_final / num_snps)
}

# Compute final SNP × ENV GRM
K_snp_env_final <- compute_snp_grm_and_add(snp_matrix_filtered, K_env, K_snp_env_final)

# Display result
head(K_snp_env_final)
