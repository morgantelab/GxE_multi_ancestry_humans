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
  make_option(c("-p", "--pcs"), type = "character", default = NULL,
              help = "Path to the scaled pcs dataset RData file", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output directory for saving results", metavar = "character"),
  make_option(c("-s", "--scratch"), type = "character", default = "/scratch3/kgoda/ukbiobank_files/tmp/snakemake_runs",
              help = "Temporary directory for storing model files", metavar = "character"),
  make_option(c("-a", "--ethn"), type = "character", default = NULL,
              help = "taking in ethn file name to get ids", metavar = "character"),
  make_option(c("-b", "--bp"), type = "character", default = NULL,
              help = "taking in bp to get type", metavar = "character"),
  make_option(c("-n", "--ethnnum"), type = "character", default = NULL,
              help = "taking in ancestry to get ancestry", metavar = "character"),
  make_option(c("-w", "--gen_eigen"), type = "character", default = NULL,
              help = "taking in pca for pcrelate grm to create W", metavar = "character"),
  make_option(c("-e", "--emat"), type = "character", default = NULL,
              help = "taking in Emat", metavar = "character")
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

# Extract `type` and `grm_source` from output file path
#output_file <- basename(opt$result)  # Get the file name from the path
#parsed_info <- strsplit(output_file, "_")[[1]]  # Split the file name by underscore

# Assuming the file name follows the pattern: "VCEm_<type>_<grm_source>_S_A_X2.csv"
#type <- parsed_info[2]  # Extract the type (e.g., sp, dp, pp)
#ethn_number <- parsed_info[4]  # Extract the ethn component without ".csv"

type <- opt$bp
ethn_number <- opt$ethnnum

# Load dataset
load(opt$data)

# Load ethn
ethn <- read.table(opt$ethn)
ethn_ids <- ethn$V1

# Load eigen results
eigen_results <- readRDS(opt$gen_eigen)
eigenvectors <- eigen_results$vectors
eigenvalues <- eigen_results$values

# Filter and scale eigenvectors by positive eigenvalues
positive_indices <- which(eigenvalues > 0)
filtered_eigenvectors <- eigenvectors[, positive_indices]
filtered_eigenvalues <- eigenvalues[positive_indices]
for(i in 1:ncol(filtered_eigenvectors)) {
  filtered_eigenvectors[, i] <- filtered_eigenvectors[, i] * sqrt(filtered_eigenvalues[i])
}
W <- filtered_eigenvectors

individual_ids <- rownames(W)
matched_dataset <- dataset[match(individual_ids, dataset$ID), ]

# Prepare phenotype vector based on type
y <- matched_dataset[[paste0(type, "0s")]]
rownames(y) <- rownames(W)

# Set the IDs in ethn to NA in y
if (!is.null(rownames(y))) {
  y[rownames(y) %in% ethn_ids] <- NA
} else {
  stop("Row names for y are missing. Ensure y has rownames corresponding to IDs.")
}

# Load scaled PCs and match to individual IDs
pcs_scaled <- readRDS(opt$pcs)
P <- pcs_scaled[match(individual_ids, pcs_scaled$ID), 2:11]

# Load eigen of E
E_eigen <- readRDS(opt$emat)
E_eigenvectors <- E_eigen$vectors
E_eigenvalues <- E_eigen$values

# Filter and scale eigenvectors by positive eigenvalues
E_positive_indices <- which(E_eigenvalues > 0)
E_filtered_eigenvectors <- E_eigenvectors[, E_positive_indices]
E_filtered_eigenvalues <- E_eigenvalues[E_positive_indices]
for(i in 1:ncol(E_filtered_eigenvectors)) {
  E_filtered_eigenvectors[, i] <- E_filtered_eigenvectors[, i] * sqrt(E_filtered_eigenvalues[i])
}
E <- E_filtered_eigenvectors

# Model setup
iter <- 90000
burnin <- 40000
thin <- 50
verb <- T
nrow_varabs <- (iter-burnin)/thin

ETA <- list(X2=list(X=P, model="FIXED", saveEffects=TRUE),
            G=list(X=W, model="BRR", saveEffects=TRUE),
            E=list(X=E, model="BRR", saveEffects=TRUE)
)

if (!is.numeric(ETA$X2$X)) {
  ETA$X2$X <- as.matrix(ETA$X2$X)
  ETA$X2$X <- apply(ETA$X2$X, 2, as.numeric)
}
print("Model ETA created")

# Run BGLR model

model <- BGLR(y=y, ETA=ETA, nIter=iter, burnIn=burnin, thin=thin, verbose=verb, saveAt=paste(opt$scratch, '/PREDs_', type, '_ethn_', ethn_number, '_S_A_X2_G_E', sep=''))

# Combine predictions and observed values into a data frame
preds <- data.frame(ID=rownames(y), Observed=model$y, Predicted=model$yHat)
write.csv(preds, file=file.path(opt$output, paste0("PREDs_", type, "_ethn_", ethn_number, "_S_A_X2_G_E.csv")), row.names=FALSE)

