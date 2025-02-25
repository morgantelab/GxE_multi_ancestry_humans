# V5 Interaction GRM ##

rm(list=ls()); gc()
set.seed(1123)

# Load required libraries
library(matrixcalc)
library(genio)
library(data.table)
library(Matrix)

# Load PLINK data using genio
plink_data <- read_plink("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/merged_geno_common_snps_selected_indiv")  # Reads data.bed, data.bim, data.fam
snp_matrix <- plink_data$X  # Extract genotype matrix
snp_subset <- snp_matrix[1:5, 1:5]  # First 5 SNPs, first 5 individuals

# Load environment data
load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/Emat_20250106.RData")
env_matrix <- Emat
env_subset <- env_matrix[1:5, 1:2]  # First 5 individuals, first 2 ENV variables

# Function to compute SNP Relationship Matrix with missing SNP handling
compute_snp_grm <- function(snp_vector) {
  n <- length(snp_vector)
  
  # Compute allele frequency while ignoring NAs
  p <- mean(snp_vector, na.rm = TRUE) / 2
  
  # Center genotype values while ignoring NAs
  centered <- snp_vector - (2 * p)
  
  # Replace NAs with the mean-centered value (imputation)
  centered[is.na(centered)] <- 0  # Conservative approach
  
  # Compute relationship matrix
  G <- outer(centered, centered) / (2 * p * (1 - p))
  return(G)
}

compute_env_grm <- function(env_vector) {
  n <- length(env_vector)
  
  # Directly compute the outer product (no need to center/scale)
  E <- outer(env_vector, env_vector) / n
  return(E)
}


# Compute SNP × ENV Interaction Matrix for selected SNPs and ENV variables
compute_snp_env_interaction <- function(snp_matrix, env_matrix) {
  num_snps <- nrow(snp_matrix)  # SNPs are now rows
  num_envs <- ncol(env_matrix)
  n <- ncol(snp_matrix)  # Individuals are now columns
  
  interaction_matrices <- list()
  
  for (s in 1:num_snps) {
    for (e in 1:num_envs) {
      # Compute SNP and ENV relationship matrices
      K_snp <- compute_snp_grm(snp_matrix[s, ])
      K_env <- compute_env_grm(env_matrix[, e])
      
      # Compute Hadamard product
      K_snp_env <- hadamard.prod(K_snp,K_env)

      # Store result
      interaction_matrices[[length(interaction_matrices) + 1]] <- K_snp_env  
    }
  }
  
  # Average over all SNP × ENV interaction matrices
  K_final <- Reduce("+", interaction_matrices) / length(interaction_matrices)
  
  return(K_final)
}

# Compute SNP × ENV Genomic Relationship Matrix with selected 5 SNPs and 2 ENV variables
K_snp_env <- compute_snp_env_interaction(snp_subset, env_subset)

# Display results
print(K_snp_env)

