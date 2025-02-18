#Interaction GRM (One Fold, All Environments) - Dynamic BP Handling ##

rm(list=ls()); gc()
set.seed(1123)

# Load necessary libraries
library(matrixcalc)
library(genio)
library(data.table)
library(Matrix)

# Define list of known interaction terms (environment terms)
interaction_terms <- c("Townsend", "act0_d", "TVtime", "sleep_d", "smoking_now", "veg_cook", "fish_oily",
                       "fish_lean", "meat_proc", "poultry", "beef", "lamb", "pork", "cheese", "salt", "tea",
                       "alc1", "waist", "getup", "coffee", "smoked_past", "BFP", "sleep_dev")

# Define parameters
fold <- 1  # Change this for other folds (1,2,3,4,5)
bp_type <- "SP"  # Options: "SP", "DP", "PP" (Dynamic BP handling)

# Define directories
plink_results_dir <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/gwas_snp_env/"
geno_data_path <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/merged_geno_common_snps_selected_indiv"
env_data_path <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/Emat_20250106.RData"

# Load PLINK data
plink_data <- read_plink(geno_data_path)
snp_matrix <- plink_data$X  # Extract genotype matrix
snp_names <- plink_data$bim$id  # SNP names

# Load environmental data
load(env_data_path)
env_matrix <- Emat

# Get all filtered files for this fold and BP type
filtered_files <- list.files(path = plink_results_dir, pattern = paste0("_", fold, ".", bp_type, "0s_filtered.txt$"), full.names = TRUE)

# Initialize final GRM for this fold and BP
num_individuals <- ncol(snp_matrix)
K_fold_grm <- matrix(0, num_individuals, num_individuals)  # Pre-allocate memory
num_envs_used <- 0  # Counter for environments used

# Loop through each known environment term
for (env_term in interaction_terms) {
  
  print(paste("🔹 Processing environment:", env_term, "for Fold:", fold, "BP:", bp_type))
  
  # Identify the corresponding filtered file for this fold, BP, and environment
  filtered_file <- filtered_files[grepl(paste0("gxe_", env_term, "_", fold, ".", bp_type, "0s_filtered.txt$"), filtered_files)]
  
  if (length(filtered_file) == 0) {
    print(paste("⚠️ No filtered SNP file found for:", env_term, "in Fold", fold, "BP:", bp_type))
    next
  }
  
  # Load SNPs from the filtered file
  gxe_results <- fread(filtered_file)
  env_snps <- gxe_results$ID  # Extract SNP column
  
  # Filter SNP matrix to keep only relevant SNPs
  snp_indices <- which(snp_names %in% env_snps)
  if (length(snp_indices) == 0) {
    print(paste("⚠️ No SNPs found in genotype data for:", env_term, "in Fold", fold, "BP:", bp_type))
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
  
  # Accumulate into the full fold GRM
  K_fold_grm <- K_fold_grm + K_snp_env
  num_envs_used <- num_envs_used + 1
}

# Compute the final full averaged GRM for this fold and BP
if (num_envs_used > 0) {
  K_fold_grm <- K_fold_grm / num_envs_used  # ✅ Averaging over environments
  fold_output_file <- file.path(plink_results_dir, paste0("Hadamard_GRM_Fold", fold, "_", bp_type, ".RData"))
  save(K_fold_grm, file = fold_output_file)
  print(paste("💾 Saved Hadamard GRM for Fold", fold, "BP:", bp_type, "to", fold_output_file))
} else {
  print(paste("⚠️ No environments were successfully processed for Fold", fold, "BP:", bp_type, ". Full fold GRM was not computed."))
}

print("✅ Processing complete for Fold", fold, "BP:", bp_type)
