# Standalone R script for running predictions without Snakemake

rm(list=ls()); gc()

# Load required libraries
library(data.table)
library(Matrix)
library(BGLR)
library(readr)
library(genio)
library(dplyr)
library(matrixcalc)
library(optparse)


# Define file paths and parameters
working_dir <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model"
data_file <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/scaled_dataset_20250106.Rdata"
output_dir <- working_dir
scratch_dir <- "/scratch3/kgoda/ukbiobank_files/tmp/snakemake_runs"
ethn_file <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/asian_ids.txt"
pcs_file <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/scaled_pcs_plink.rds"
eigen_file <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/pca_for_pcrelate.rds"
envvars_file <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/E_eigen_G_E.rds"
GEselect_file <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/gwas_snp_env/Hadamard_GRM_Ancestryasian_SP_0.001.RData"

# Output file name
output_file <- paste0(output_dir, "/PREDs_SP_ethn_asianpval0.001_X1_X2_G_E_GEselect_check.csv")

# Set working directory
setwd(working_dir)

# Load dataset
load(data_file)

# Load ethnicity data
ethn <- read.table(ethn_file)
ethn_ids <- ethn$V1  # Extract only the first column

# Load eigen results
eigen_results <- readRDS(eigen_file)
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

# Prepare phenotype vector
y <- matched_dataset$SP0s
rownames(y) <- rownames(W)
# Set the IDs in ethn_1 to NA in y
if (!is.null(rownames(y))) {
  y[rownames(y) %in% ethn_ids] <- NA
} else {
  stop("Row names for y are missing. Ensure y has rownames corresponding to IDs.")
}

# Extract fixed covariates (X)
X <- dataset[, c("AOPs", "AOPss", "Sex_SI")]
X$AOPs <- as.vector(X$AOPs)
X$AOPss <- as.vector(X$AOPss)
rownames(X) <- rownames(W)

print("initial X created")

# Load and match scaled PCs
pcs_scaled <- readRDS(pcs_file)
P <- pcs_scaled[match(individual_ids, pcs_scaled$ID), 2:11]

# Load and process environmental eigenvectors
E_eigen <- readRDS(envvars_file)
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

# Load and process GEselect eigenvectors
load(GEselect_file)
GEselect_eigen <- eigen(K_ancestry_grm)
rownames(GEselect_eigen$vectors) <- rownames(K_ancestry_grm)

GEselect_eigenvectors <- GEselect_eigen$vectors
GEselect_eigenvalues <- GEselect_eigen$values

# Filter and scale eigenvectors by positive eigenvalues
GEselect_positive_indices <- which(GEselect_eigenvalues > 0)
GEselect_filtered_eigenvectors <- GEselect_eigenvectors[, GEselect_positive_indices]
GEselect_filtered_eigenvalues <- GEselect_eigenvalues[GEselect_positive_indices]
for(i in 1:ncol(GEselect_filtered_eigenvectors)) {
  GEselect_filtered_eigenvectors[, i] <- GEselect_filtered_eigenvectors[, i] * sqrt(GEselect_filtered_eigenvalues[i])
}
GEselect <- GEselect_filtered_eigenvectors

# Model parameters
iter <- 90000
burnin <- 40000
thin <- 50
verb <- T
nrow_varabs <- (iter-burnin)/thin
scratch <- "/scratch3/kgoda/ukbiobank_files/tmp/snakemake_runs" # Use this!

# Define model
ETA <- list(X1=list(X=X, model="FIXED", saveEffects=TRUE),
            X2=list(X=P, model="FIXED", saveEffects=TRUE),
            G=list(X=W, model="BRR", saveEffects=TRUE),
            E=list(X=E, model="BRR", saveEffects=TRUE),
            GxEselect=list(X=GEselect, model="BRR", saveEffects=TRUE)
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

# Run BGLR model
model <- BGLR(y = y, ETA = ETA, nIter = iter, burnIn = burnin, thin = thin, verbose = verb, 
              saveAt = paste0(scratch_dir, "/PREDs_SP_ethn_asian_X1_X2_G_E_GEselect"))

# Save predictions
preds <- data.frame(ID = rownames(y), Observed = model$y, Predicted = model$yHat)
write.csv(preds, file = output_file, row.names = FALSE)
