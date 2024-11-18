rm(list=ls()); gc()

# Load required libraries
library(data.table)
library(Matrix)
library(BGLR)
library(readr)
library(genio)
library(dplyr)

setwd("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model")

load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/scaled_dataset_20241025.Rdata")

eigen_path_pcrelate <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/pca_for_pcrelate.rds"
eigen_results_pcrelate <- readRDS(eigen_path_pcrelate)
eigenvectors <- eigen_results_pcrelate$vectors
eigenvalues <- eigen_results_pcrelate$values

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

### Extract the phenotype vectors ###
y_PP_pcrelate <- matched_dataset$PP0s

# Prepare covariate matrix
X <- matched_dataset[, c("AOPs", "AOPss", "Sex_SI")]
X$AOPs <- as.vector(X$AOPs)
X$AOPss <- as.vector(X$AOPss)
rownames(X) <- rownames(W)

# Load scaled PCs and match to individual IDs
pcs_scaled <- readRDS("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/scaled_pcs_plink.rds")
matched_pcs <- pcs_scaled[match(individual_ids, pcs_scaled$ID), 2:11]
P <- matched_pcs

print("Combined X matrix created")

# Model setup
iter <- 90000
burnin <- 40000
thin <- 50
verb <- T
nrow_varabs <- (iter-burnin)/thin
scratch <- "/scratch3/kgoda/ukbiobank_files/tmp/snakemake_runs" # Use this!

ETA <- list(X1=list(X=X, model="FIXED", saveEffects=TRUE), X2=list(X=P, model="FIXED", saveEffects=TRUE), G=list(X=W, model="BRR", saveEffects=TRUE))

if (!is.numeric(ETA$X1$X)) {
  ETA$X1$X <- as.matrix(ETA$X1$X)
  ETA$X1$X <- apply(ETA$X1$X, 2, as.numeric)
}
print("Model ETA created")

if (!is.numeric(ETA$X2$X)) {
  ETA$X2$X <- as.matrix(ETA$X2$X)
  ETA$X2$X <- apply(ETA$X2$X, 2, as.numeric)
}
print("Model ETA created")

# Run BGLR model

model <- BGLR(y=y_PP_pcrelate, ETA=ETA, nIter=iter, burnIn=burnin, thin=thin, verbose=verb, saveAt=paste(scratch, '/PP_pcrelate_pcs_plink_G_run_X_separated_', sep=''))


### Collect results ###

### Model parameters ###  

### intercept ###
zz0 <- read.table(paste(scratch, '/PP_pcrelate_pcs_plink_G_run_X_separated_mu.dat', sep=''), header=F); colnames(zz0) <- "int"

print("zz0 created")

### variance of ridge regression betas for the additive genetic effect ###
zz1 <- read.table(paste(scratch, '/PP_pcrelate_pcs_plink_G_run_X_separated_ETA_G_varB.dat', sep=''), header=F); colnames(zz1) <- "G"

print("zz1 created")

### variance of the residual deviations ###
zz9 <- read.table(paste(scratch, '/PP_pcrelate_pcs_plink_G_run_X_separated_varE.dat', sep=''), header=F); colnames(zz9) <- "res"

print("zz9 created")

### dataframe with model parameters ###
VCEm_PP_pcrelate <- data.frame(zz0, zz1, zz9)

print("PP_pcrelate_pcs_plink_G_run_X_separated created now saving")

write.csv(VCEm_PP_pcrelate, file="/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/VCEm_PP_pcrelate_pcs_plink_G_run_X_separated.csv", row.names=TRUE)
#VCEm_PP <- read.csv("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/apr_run_with_chinese_20240402/5_variance_components/VCEm_PP_pcrelate_pcs_plink_G_run_X_separated")
print("PP_pcrelate_pcs_plink_G_run_X_separated written")

### Sampled regression effects ###  

### intercept and fixed effects part 1 ###
B1 <- read.table(paste(scratch, '/PP_pcrelate_pcs_plink_G_run_X_separated_ETA_X1_b.dat', sep=''), header=T)

print("B1 done now B2")

### additive genetic random effects ###
B2 <- readBinMat(paste(scratch, '/PP_pcrelate_pcs_plink_G_run_X_separated_ETA_G_b.bin', sep=''))

### fixed effects part 2 pcs ###
#B3 <- readBinMat(paste(scratch, '/PP_pcrelate_pcs_plink_G_run_X_separated_ETA_X2_b.dat', sep=''))
B3 <- as.matrix(fread(paste(scratch, '/PP_pcrelate_pcs_plink_G_run_X_separated_ETA_X2_b.dat', sep='')))

### dataframe with variance partition ###
varabs <- matrix(NA, nrow_varabs, 3); colnames(varabs) <- c("V_X1", "V_X2", "V_G")

print("filling up cols of varab")

### fill up column for fixed covariates ###
varabs[, 1] <- matrixStats::colVars(ETA$X1$X%*%t(B1))[-c(1:(burnin/thin))]   
print("varab col one filled up now col 2")

varabs[, 2] <- matrixStats::colVars(ETA$X2$X%*%t(B3))[-c(1:(burnin/thin))]

### fill up column for random additive genetic covariates ###
varabs[, 3] <- matrixStats::colVars(tcrossprod(ETA$G$X,B2))

print("varab cols done now saving varab")

write.csv(varabs, file="/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_PP_pcrelate_pcs_plink_G_run_X_separated.csv", row.names=TRUE)

print("varab saved")

### Note -c(1:(burnin/thin)). This is because burn-in samples are stored in B1, but are not in B2 and B3.

