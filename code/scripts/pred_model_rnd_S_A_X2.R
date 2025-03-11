rm(list=ls()); gc()
set.seed(1123)

# Load required libraries
library(data.table)
library(Matrix)
library(BGLR)
library(readr)
library(genio)
library(dplyr)
library(optparse)

# Command-line options
option_list <- list(
  make_option(c("-d", "--dir"), type = "character", default = NULL,
              help = "path to the working directory", metavar = "character"),
  make_option(c("-t", "--data"), type = "character", default = NULL,
              help = "Path to the scaled dataset RData file", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output directory for saving results", metavar = "character"),
  make_option(c("-s", "--scratch"), type = "character", default = "/scratch3/kgoda/ukbiobank_files/tmp/snakemake_runs",
              help = "Temporary directory for storing model files", metavar = "character"),
  make_option(c("-f", "--fold"), type = "character", default = NULL,
              help = "taking in fold file name to get ids", metavar = "character"),
  make_option(c("-p", "--pcs"), type = "character", default = NULL,
              help = "taking in pcs file name to get scaled pcs", metavar = "character"),
  make_option(c("-b", "--bp"), type = "character", default = NULL,
              help = "taking in bp to get type", metavar = "character"),
  make_option(c("-n", "--num"), type = "character", default = NULL,
              help = "taking in num to get fold number", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check required arguments
if (is.null(opt$data) | is.null(opt$pcs) | is.null(opt$output)) {
  print_help(opt_parser)
  stop("All file paths and phenotype type must be provided.", call. = FALSE)
}

# Set the working directory dynamically
setwd(opt$dir)

# # Parse `type` and `fold` from the result file name
# output_file <- basename(opt$result)  # Get the file name from the path
# parsed_info <- strsplit(output_file, "_")[[1]]  # Split the file name by underscore
#
# type <- parsed_info[2]  # Extract the type component
# fold_number <- parsed_info[4]  # Extract the fold component without ".csv"

type <- opt$bp
fold_number <- opt$num

# Load dataset
load(opt$data)

# Load fold
fold <- readRDS(opt$fold)
fold_ids <- fold$ID

# Prepare phenotype vector based on type
### Extract the phenotype vectors ###
y <- dataset[[paste0(type, "0s")]]
rownames(y) <- dataset$ID

# Set the IDs in fold_1 to NA in y
if (!is.null(rownames(y))) {
  y[rownames(y) %in% fold_ids] <- NA
} else {
  stop("Row names for y are missing. Ensure y has rownames corresponding to IDs.")
}

# Load scaled PCs and match to individual IDs
pcs_scaled <- readRDS(opt$pcs)
P <- pcs_scaled[match(dataset$ID, pcs_scaled$ID), 2:11]

# Model setup
iter <- 90000
burnin <- 40000
thin <- 50
verb <- T
nrow_varabs <- (iter-burnin)/thin

ETA <- list(X2=list(X=P, model="FIXED", saveEffects=TRUE))

if (!is.numeric(ETA$X2$X)) {
  ETA$X2$X <- as.matrix(ETA$X2$X)
  ETA$X2$X <- apply(ETA$X2$X, 2, as.numeric)
}
print("Model ETA created")

# Run BGLR model

model <- BGLR(y=y, ETA=ETA, nIter=iter, burnIn=burnin, thin=thin, verbose=verb, saveAt=paste(opt$scratch, '/PREDs_', type, '_fold_', fold_number, '_S_A_X2_', sep=''))

# Combine predictions and observed values into a data frame
preds <- data.frame(ID=rownames(y), Observed=model$y, Predicted=model$yHat)
write.csv(preds, file=file.path(opt$output, paste0("PREDs_", type, "_Fold_", fold_number, "_S_A_X2.csv")), row.names=FALSE)


