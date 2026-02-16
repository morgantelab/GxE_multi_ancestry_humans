# Clear environment
rm(list = ls())

set.seed(1123)

# Load required libraries
library(dplyr)
library(optparse)

# Define command line arguments
option_list <- list(
  make_option(c("--input_rds"), type = "character", help = "Path to RDS file with dataset"),
  make_option(c("--output_eigen"), type = "character", help = "Output RDS file for E_eigen"),
  make_option(c("--output_Emat"), type = "character", help = "Output RData file for scaled E matrix")
)

# Parse arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Load dataset
dataset <- readRDS(opt$input_rds)

# Define environment variables
Evars <- c("Townsend", "act0_d", "TVtime", "sleep_d", "smoking_now", "veg_cook", "fish_oily",
           "fish_lean", "meat_proc", "poultry", "beef", "lamb", "pork", "cheese", "salt", "tea",
           "alc1", "waist", "getup", "coffee", "smoked_past", "BFP", "sleep_dev")

# Convert to matrix and scale
Emat <- as.matrix(dataset[Evars])
rownames(Emat) <- dataset$ID
Emat <- scale(Emat)

# Compute environmental kernel
E <- tcrossprod(Emat)
E <- E / mean(diag(E))

# Compute eigen-decomposition
E_eigen <- eigen(E)
rownames(E_eigen$vectors) <- rownames(Emat)

# Save outputs
saveRDS(E_eigen, file = opt$output_eigen)
save(Emat, file = opt$output_Emat)
