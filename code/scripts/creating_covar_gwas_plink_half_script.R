rm(list=ls()); gc()
library(data.table)
library(Matrix)

# Load dataset for phenotype
load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/scaled_dataset_20250106.Rdata")
# Extract required columns: ID, DP0s, SP0s, PP0s
phenotype <- data.table(ID = as.character(dataset$ID), age = as.numeric(dataset$AOPs), age2 = as.numeric(dataset$AOPss), sex = as.factor(dataset$Sex_SI), DP0s = as.numeric(dataset$DP0s), SP0s = as.numeric(dataset$SP0s), PP0s = as.numeric(dataset$PP0s))

# Load principal components (PCs) from PLINK
pcs_scaled <- readRDS("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/scaled_pcs_plink.rds")
# Convert `ID` in `pcs_scaled` to character
pcs_scaled$ID <- as.character(pcs_scaled$ID)

# Merge phenotype, PCs
merged_pheno_pcs <- Reduce(function(x, y) merge(x, y, by="ID", all=TRUE), list(phenotype, pcs_scaled))

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

# Define output file path
output_file <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/gwas_snp_env/covar_pheno.txt"

# Ensure all required columns exist and are properly formatted
plink_data <- final_data[, .(ID, ID, DP0s, SP0s, PP0s, age, age2, sex, X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, Townsend, act0_d, TVtime, sleep_d, smoking_now, veg_cook, fish_oily, fish_lean, meat_proc, poultry, beef, lamb, pork, cheese, salt, tea, alc1, waist, getup, coffee, smoked_past, BFP, sleep_dev)]

# Rename first two columns to match PLINK's required format (FID, IID)
names(plink_data)[1:2] <- c("FID", "IID")

# Save as tab-separated file
fwrite(plink_data, output_file, sep="\t", quote=FALSE, col.names=TRUE)

# Confirm file is written
print(paste("File saved:", output_file))
