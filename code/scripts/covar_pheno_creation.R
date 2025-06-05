rm(list=ls()); gc()
set.seed(1123)

# Load libraries
library(data.table)
library(Matrix)
library(optparse)

# Define command line options
option_list <- list(
  make_option(c("--data3"), type = "character", help = "Path to data3 RDS file"),
  make_option(c("--pcs"), type = "character", help = "Path to scaled PLINK PCs RDS file"),
  make_option(c("--emat"), type = "character", help = "Path to Emat RData file"),
  make_option(c("--output"), type = "character", help = "Path to output covar_pheno.txt")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Load data3 RDS
dataset <- readRDS(opt$data3)
phenotype <- data.table(ID = as.character(dataset$ID),
                        age = as.numeric(dataset$AOPs),
                        age2 = as.numeric(dataset$AOPsss),
                        sex = as.numeric(dataset$Sex_SIs),
                        DP0s = as.numeric(dataset$DP0s),
                        SP0s = as.numeric(dataset$SP0s),
                        PP0s = as.numeric(dataset$PP0s))

# Load PCs
pcs_scaled <- readRDS(opt$pcs)
pcs_scaled$ID <- as.character(pcs_scaled$ID)

# Merge phenotype and PCs
merged_pheno_pcs <- merge(phenotype, pcs_scaled, by="ID", all=TRUE)

# Load Emat
load(opt$emat)
if (!exists("Emat")) stop("Emat object not found in loaded RData")

# Convert Emat to data.table
env <- as.data.table(Emat, keep.rownames="ID")
env[, ID := as.character(ID)]

# Merge environment with phenotype+PCs
final_data <- merge(merged_pheno_pcs, env, by="ID", all=TRUE)
if (nrow(final_data) != nrow(merged_pheno_pcs)) stop("Mismatch after merging environment data")

# Format output
plink_data <- final_data[, .(FID = ID, IID = ID, DP0s, SP0s, PP0s, age, age2, sex, X1, X2, X3, X4, X5, X6, X7, X8, X9, X10,
                             Townsend, act0_d, TVtime, sleep_d, smoking_now, veg_cook, fish_oily, fish_lean,
                             meat_proc, poultry, beef, lamb, pork, cheese, salt, tea, alc1, waist, getup, coffee,
                             smoked_past, BFP, sleep_dev)]

# Write output
fwrite(plink_data, opt$output, sep="\t", quote=FALSE, col.names=TRUE)
cat("File saved:", opt$output, "\n")

