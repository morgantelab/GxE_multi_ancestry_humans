rm(list=ls()); gc()
library(data.table)
library(Matrix)
library(lmtest)
library(qqman)
library(genio)

# Load precomputed GxE interaction matrix (PC-Relate)
#load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/GE_hadamard_prod_pcrelate.RData")  # Loads GE_pcrelate
# Extract individual-level GxE interaction effects (diagonal of GE_pcrelate)
#diagGE <- diag(GE_pcrelate)
#diagGE <- data.table(ID = rownames(GE_pcrelate), diagGE = diagGE)

# Load dataset for phenotype
load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/scaled_dataset_20250106.Rdata")
# Extract required columns: ID, DP0s, SP0s, PP0s
phenotype <- data.table(ID = as.character(dataset$ID), age = as.numeric(dataset$AOPs), age2 = as.numeric(dataset$AOPss), sex = as.numeric(dataset$Sex_SI), DP0s = as.numeric(dataset$DP0s), SP0s = as.numeric(dataset$SP0s), PP0s = as.numeric(dataset$PP0s))

# Load principal components (PCs) from PLINK
pcs_scaled <- readRDS("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/scaled_pcs_plink.rds")
# Convert `ID` in `pcs_scaled` to character
pcs_scaled$ID <- as.character(pcs_scaled$ID)

# Merge phenotype, PCs
merged_pheno_pcs <- Reduce(function(x, y) merge(x, y, by="ID", all=TRUE), list(phenotype, pcs_scaled))

# Ensure GE_pcrelate has the same individuals as `data`
#GE_pcrelate <- GE_pcrelate[match(data$ID, rownames(GE_pcrelate)), , drop=FALSE]
#if (nrow(GE_pcrelate) != nrow(data)) stop("Error: Mismatch in sample IDs between GE_pcrelate and phenotype data.")

# Extract PC column names
pc_names <- names(pcs_scaled)[-1]  # Exclude ID column

# Load environment data
load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/Emat_20250106.RData")

# Convert environment data into a data.table
env <- as.data.table(Emat, keep.rownames="ID")  # Move row names to `ID` column
env[, ID := as.character(ID)]  # Ensure ID is character

# Merge environment data with phenotype & PCs
final_data <- merge(merged_pheno_pcs, env, by="ID", all=TRUE)

# Ensure final dataset has the correct number of individuals
if (nrow(final_data) != nrow(merged_pheno_pcs)) stop("Error: Mismatch after merging environment data.")


# Load SNP genotype data from PLINK (.bed, .bim, .fam)
plink_snp <- read_plink("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/merged_geno_common_snps_selected_indiv")
#bed_file <- read_plink("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/merged_geno_common_snps_selected_indiv.bed")
#bim_file <- read_plink("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/merged_geno_common_snps_selected_indiv.bim")
#fam_file <- read_plink("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/merged_geno_common_snps_selected_indiv.fam")

#working on rest of this. 30th Jan 2025

# Extract genotype matrix
geno_matrix <- plink_snp$X  # Matrix of SNPs (rows: individuals, columns: SNPs)
geno_matrix <- t(geno_matrix)  # Transpose: now individuals are rows
geno_ids <- plink_snp$fam$id  # Extract individual IDs
snp_names <- plink_snp$bim$id # Extract snp ids

# Ensure SNP matrix matches the phenotype dataset
#geno_matrix <- geno_matrix[match(merged_pheno_pcs$ID, geno_ids), , drop=FALSE]
#if (nrow(geno_matrix) != nrow(merged_pheno_pcs)) stop("Error: Mismatch between SNP genotype data and phenotype data.")

# Select environment for testing
env_name <- "Townsend"
trait <- "SP0s"  # Change to test different traits

# Assign environment data
final_data[, env := get(env_name)]  # Assign the correct environment column

# Function to run GWAS for each SNP
run_gxe_gwas <- function(snp_idx) {
  snp_data <- geno_matrix[ , snp_idx]  # Extract SNP values
  final_data[, snp := snp_data]  # Assign SNP genotype to dataset

  # Define model formula using automatic interaction syntax
  formula <- as.formula(paste(trait, "~ age + I(age^2) + sex +",
                              paste(pc_names, collapse=" + "),
                              "+ snp * env"))  # SNP x Env interaction

  # Fit linear model
  model <- lm(formula, data=final_data, na.action=na.omit)

  # Extract coefficients for SNP × Env interaction
  coef_summary <- summary(model)$coefficients

  # Store results for SNP × Env interaction
  return(data.table(
    Trait = trait,
    SNP = snp_names[snp_idx],
    Environment = env_name,
    Beta = coef_summary["snp:env", "Estimate"],
    SE = coef_summary["snp:env", "Std. Error"],
    P = coef_summary["snp:env", "Pr(>|t|)"]
  ))
}

# Run GWAS for all SNPs using lapply for efficiency
results_list <- lapply(seq_len(ncol(geno_matrix)), run_gxe_gwas)

# Combine results into a single data.table
results <- rbindlist(results_list)

# Save results
fwrite(results, "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/gwas_snp_env/gxe_SP0s_townsend_all_snps_lm_results.txt", sep="\t")

# Print summary of results
print(head(results))





#
# # Select a single SNP and a single environment for testing
# snp_idx <- 1  # Change this to test different SNPs
# env_name <- "Townsend"  # Change this to test different environments
# trait <- "DP0s"  # Change this to test different traits
#
# # Get SNP genotype data
# snp_data <- geno_matrix[, snp_idx]
#
# # Add SNP and environment data to dataset
# final_data$snp <- snp_data
# final_data$env <- final_data[[env_name]]
#
# # Define model formula using automatic interaction syntax
# formula <- as.formula(paste(trait, "~ age + I(age^2) + sex +",
#                             paste(pc_names, collapse=" + "),
#                             "+ snp * env"))  # SNP x Env interaction handled automatically
#
# # Fit linear model
# model <- lm(formula, data=final_data, na.action=na.omit)
#
# # Extract coefficients for SNP × Env interaction
# coef_summary <- summary(model)$coefficients
#
# # Store results in a data table
# results <- data.table(
#   Trait = trait,
#   SNP = snp_names[snp_idx],
#   Environment = env_name,
#   Beta = coef_summary["snp:env", "Estimate"],
#   SE = coef_summary["snp:env", "Std. Error"],
#   P = coef_summary["snp:env", "Pr(>|t|)"]
# )
#
# # Save results
# fwrite(results, "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/results/single_test_gxe_lm_results.txt", sep="\t")
#
# # Print results
# print(results)
