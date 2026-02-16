rm(list=ls()); gc()
set.seed(1123)

# Load required libraries
.libPaths('/data2/morgante_lab/kgoda/software/R/x86_64-pc-linux-gnu-library/4.2')

# Load libraries
library(data.table)
library(Matrix)
library(optparse)

# Define command line options
option_list <- list(
  make_option(c("--data"), type = "character", help = "Path to data3 RDS file"),
  make_option(c("--pcs"), type = "character", help = "Path to scaled PLINK PCs RDS file"),
  make_option(c("--emat"), type = "character", help = "Path to Emat RData file"),
  make_option(c("--output"), type = "character", help = "Path to output covar_pheno.txt")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if(0){
  # Dummy loading for testing
  load('data/white_24k/Emat_20251104.RData')
  dataset <- readRDS('data/white_24k/data4_20251104.rds')
  pcs_scaled <- readRDS('data/white_24k/scaled_pcs_plink.rds')
}

# Load data3 RDS
dataset <- readRDS(opt$data)
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

envlist <- c('Townsend', 'act0_d', 'TVtime', 'sleep_d', 'smoking_now', 'veg_cook', 'fish_oily', 'fish_lean',
             'meat_proc', 'poultry', 'beef', 'lamb', 'pork', 'cheese', 'salt', 'tea', 'alc1', 'waist', 'getup', 'coffee',
             'smoked_past', 'BFP', 'sleep_dev')

# Define the column groups
set1 <- envlist
set2 <- paste0("X",1:10)

# Straight from chat_Gepeto to merge sets of columns - nklimko
# Compute all pairwise products
product_list <- lapply(set1, function(col1) {
  lapply(set2, function(col2) {
    new_col <- plink_data[[col1]] * plink_data[[col2]]
    name <- paste(col1, col2, sep = "__")
    setNames(data.frame(new_col), name)
  }) |> do.call(what = cbind)
}) |> do.call(what = cbind)

# Append results to the existing dataframe
final_data <- cbind(plink_data, product_list)

# names(final_data)


# Write output
fwrite(final_data, opt$output, sep="\t", quote=FALSE, col.names=TRUE)
cat("File saved:", opt$output, "\n")

