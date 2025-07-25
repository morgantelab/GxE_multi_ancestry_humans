set.seed(1123)
### Running Predictions for model with just X ###
rm(list=ls()); gc()

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
  make_option(c("-r", "--result"), type = "character", default = NULL,
              help = "taking in output file name to get type", metavar = "character"),
  make_option(c("-a", "--ethn"), type = "character", default = NULL,
              help = "taking in ethn file name to get ids", metavar = "character"),
  make_option(c("-p", "--pcs"), type = "character", default = NULL,
              help = "taking in pcs file name to get scaled pcs", metavar = "character"),
  make_option(c("-e", "--eigen"), type = "character", default = NULL,
              help = "taking in eigen of pcrelate grm file name to get g", metavar = "character"),
  make_option(c("-v", "--envvars"), type = "character", default = NULL,
              help = "envvar matrix", metavar = "character"),
  make_option(c("-q", "--GE"), type = "character", default = NULL,
              help = "Path to the scaled PCs RDS file", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Set the working directory dynamically
setwd(opt$dir)

# Parse `type` and `ethn` from the result file name
output_file <- basename(opt$result)  # Get the file name from the path
parsed_info <- strsplit(output_file, "_")[[1]]  # Split the file name by underscore

type <- parsed_info[2]  # Extract the type component
ethn_number <- parsed_info[4]  # Extract the ethn component without ".csv"

# Load dataset
dataset <- readRDS(opt$data)

# Load ethn
ethn <- read.table(opt$ethn)
ethn_ids <- ethn$V1

# Load eigen results
eigen_results <- readRDS(opt$eigen)
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

# Set the IDs in ethn_1 to NA in y
if (!is.null(rownames(y))) {
  y[rownames(y) %in% ethn_ids] <- NA
} else {
  stop("Row names for y are missing. Ensure y has rownames corresponding to IDs.")
}

### X is the incidence matrix for the 'fixed' covariates (no penalisation, no shrinkage). here age and sex ###
### Extract the covariate matrix ###
X <- dataset[, c("AOPs", "AOPsss", "Sex_SIs")]
X$AOPs <- as.vector(X$AOPs)
X$AOPsss <- as.vector(X$AOPsss)
X$Sex_SIs <- as.vector(X$Sex_SIs)
rownames(X) <- rownames(W)

print("initial X created")

# Load scaled PCs and match to individual IDs
pcs_scaled <- readRDS(opt$pcs)
matched_pcs <- pcs_scaled[match(individual_ids, pcs_scaled$ID), 2:11]
P <- matched_pcs

# Load eigen of E
E_eigen <- readRDS(opt$envvars)
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

# Load eigen of GE
GE_eigen <- readRDS(opt$GE)
GE_eigenvectors <- GE_eigen$vectors
GE_eigenvalues <- GE_eigen$values

# Filter and scale eigenvectors by positive eigenvalues
GE_positive_indices <- which(GE_eigenvalues > 0)
GE_filtered_eigenvectors <- GE_eigenvectors[, GE_positive_indices]
GE_filtered_eigenvalues <- GE_eigenvalues[GE_positive_indices]
for(i in 1:ncol(GE_filtered_eigenvectors)) {
  GE_filtered_eigenvectors[, i] <- GE_filtered_eigenvectors[, i] * sqrt(GE_filtered_eigenvalues[i])
}
GE <- GE_filtered_eigenvectors

### Model ###

iter <- 90000
burnin <- 40000
thin <- 50
verb <- T
nrow_varabs <- (iter-burnin)/thin
scratch <- "/scratch3/kgoda/ukbiobank_files/tmp/snakemake_runs" # Use this!

### Define model ###

ETA <- list(X1=list(X=X, model="FIXED", saveEffects=TRUE),
            X2=list(X=P, model="FIXED", saveEffects=TRUE),
            G=list(X=W, model="BRR", saveEffects=TRUE),
            E=list(X=E, model="BRR", saveEffects=TRUE),
            GxE=list(X=GE, model="BRR", saveEffects=TRUE)
)

### Ensure ETA$X$X is numeric ###
if (!is.numeric(ETA$X1$X)) {
  ETA$X1$X <- as.matrix(ETA$X1$X)
  ETA$X1$X <- apply(ETA$X1$X, 2, as.numeric)
}

if (!is.numeric(ETA$X2$X)) {
  ETA$X2$X <- as.matrix(ETA$X2$X)
  ETA$X2$X <- apply(ETA$X2$X, 2, as.numeric)
}

print("model ETA created")

model <- BGLR(y=y, ETA=ETA, nIter=iter, burnIn=burnin, thin=thin, verbose=verb, saveAt=paste(opt$scratch, '/PREDs_', type, '_ethn_', ethn_number, '_X1_X2_G_E_GE', sep=''))

# Combine predictions and observed values into a data frame
preds <- data.frame(ID=rownames(y), Observed=model$y, Predicted=model$yHat)
write.csv(preds, file=file.path(opt$output, paste0("PREDs_", type, "_ethn_", ethn_number, "_X1_X2_G_E_GE.csv")), row.names=FALSE)
