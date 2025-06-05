### For PCRelate and plink GRM ###
rm(list=ls()); gc()

# Load required libraries
library(data.table)
library(Matrix)
library(BGLR)
library(readr)
library(genio)
library(dplyr)
library(optparse)


### Please refer to https://github.com/gdlc/BGLR-R for a more valuable description ###

setwd("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model")

### load dataset ###
load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/scaled_dataset_20241025.Rdata")
print("dataset loaded")

# Load the eigen results
eigen_path_plink <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/pca_for_plink.rds"
eigen_results_plink <- readRDS(eigen_path_plink)
eigenvectors <- eigen_results_plink$vectors
eigenvalues <- eigen_results_plink$values
print("eigen results handled")

### Filter out negative eigenvalues and those close to zero ###
positive_indices <- which(eigenvalues > 0)

### Filter the eigenvectors and eigenvalues ###
filtered_eigenvectors <- eigenvectors[, positive_indices]
filtered_eigenvalues <- eigenvalues[positive_indices]

### Scale eigenvectors by the square root of their corresponding positive eigenvalues ###
for(i in 1:ncol(filtered_eigenvectors)){
  filtered_eigenvectors[,i] <- filtered_eigenvectors[,i] * sqrt(filtered_eigenvalues[i])
}

### Assign the filtered and scaled eigenvectors to W ###
W <- filtered_eigenvectors

print("W created")

### y is the vector of phenotypes ###
### rownames for phenotype vector ###
individual_ids <- rownames(W)

### Ensure that dataset is sorted in the order of individual_ids ###
matched_dataset <- dataset[match(individual_ids, dataset$ID), ]

### Extract the phenotype vectors ###
y_PP_plink <- matched_dataset$PP0s

print("y_PP_plink created")

### X is the incidence matrix for the 'fixed' covariates (no penalisation, no shrinkage). here age and sex ###
### Extract the covariate matrix ###
X <- matched_dataset[, c("AOPs", "AOPss", "Sex_SI")]
X$AOPs <- as.vector(X$AOPs)
X$AOPss <- as.vector(X$AOPss)
rownames(X) <- rownames(W)

print("initial X created")

E_eigen_path <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/E_eigen_G_E.rds"
E_eigen_results <- readRDS(E_eigen_path)
E_eigenvectors <- E_eigen_results$vectors
E_eigenvalues <- E_eigen_results$values

# Filter and scale eigenvectors by positive eigenvalues
E_positive_indices <- which(E_eigenvalues > 0)
E_filtered_eigenvectors <- E_eigenvectors[, E_positive_indices]
E_filtered_eigenvalues <- E_eigenvalues[E_positive_indices]
for(i in 1:ncol(E_filtered_eigenvectors)) {
  E_filtered_eigenvectors[, i] <- E_filtered_eigenvectors[, i] * sqrt(E_filtered_eigenvalues[i])
}
E <- E_filtered_eigenvectors

# Load eigen of GE
GE_eigen <- readRDS("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/eigen_GE_plink.rds")
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

# Model setup
iter <- 90000
burnin <- 40000
thin <- 50
verb <- T
nrow_varabs <- (iter-burnin)/thin
scratch <- "/scratch3/kgoda/ukbiobank_files/tmp/snakemake_runs" # Use this!


ETA <- list(X1=list(X=X, model="FIXED", saveEffects=TRUE),
            G=list(X=W, model="BRR", saveEffects=TRUE),
            E=list(X=E, model="BRR", saveEffects=TRUE),
            GxE=list(X=GE, model="BRR", saveEffects=TRUE)
)

if (!is.numeric(ETA$X1$X)) {
  ETA$X1$X <- as.matrix(ETA$X1$X)
  ETA$X1$X <- apply(ETA$X1$X, 2, as.numeric)
}
print("Model ETA created")

# Run BGLR model

model <- BGLR(y=y, ETA=ETA, nIter=iter, burnIn=burnin, thin=thin, verbose=verb, saveAt=paste(opt$scratch, '/PP_plink_GE_2corr_' , sep=''))

###Collecting results###
# Load results from BGLR output
zz0 <- read.table(paste(scratch, '/PP_plink_GE_2corr_mu.dat', sep=''), header=F)
colnames(zz0) <- "int"

zz1 <- read.table(paste(scratch, '/PP_plink_GE_2corr_ETA_G_varB.dat', sep=''), header=F)
colnames(zz1) <- "G"

zz2 <- read.table(paste(scratch, '/PP_plink_GE_2corr_ETA_E_varB.dat', sep=''), header=F)
colnames(zz2) <- "E"

zz9 <- read.table(paste(scratch, '/PP_plink_GE_2corr_varE.dat', sep=''), header=F)
colnames(zz9) <- "res"

# Combine results into a dataframe and save
VCEm <- data.frame(zz0, zz1, zz2, zz9)
write.csv(VCEm, file="/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/VCEm_PP_plink_GE.csv", row.names=TRUE)

# Sampled regression effects
B1 <- read.table(paste(scratch, '/PP_plink_GE_2corr_ETA_X1_b.dat', sep=''), header=TRUE)
B2 <- readBinMat(paste(scratch, '/PP_plink_GE_2corr_ETA_G_b.bin', sep=''))
B3 <- readBinMat(paste(scratch, '/PP_plink_GE_2corr_ETA_E_b.bin', sep=''))
B4 <- readBinMat(paste(scratch, '/PP_plink_GE_2corr_ETA_GxE_b.bin', sep=''))

# Calculate variance components
varabs <- matrix(NA, nrow_varabs, 4); colnames(varabs) <- c("V_X1", "V_G", "V_E", "V_GxE")

# Fill variance components
varabs[, 1] <- matrixStats::colVars(ETA$X1$X %*% t(B1))[-c(1:(burnin/thin))]
varabs[, 2] <- matrixStats::colVars(tcrossprod(ETA$G$X, B2))
varabs[, 3] <- matrixStats::colVars(tcrossprod(ETA$E$X, B3))
varabs[, 4] <- matrixStats::colVars(tcrossprod(ETA$GxE$X, B4))

# Save variance components
write.csv(varabs, file="/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_PP_plink_GE.csv", row.names=TRUE)
print("Variance partition results saved")
