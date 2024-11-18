### For PCRelate and plink GRM ###
rm(list=ls()); gc()

library(data.table)
library(Matrix)
library(BGLR)
library(readr)
library(genio)
library(dplyr)

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

E_eigen_path <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/E_eigen_DP_plink_G_E.rds"
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

### Model ###

iter <- 90000
burnin <- 40000
thin <- 50
verb <- T
nrow_varabs <- (iter-burnin)/thin
scratch <- "/scratch3/kgoda/ukbiobank_files/tmp/snakemake_runs" # Use this!

### Define model ###

ETA <- list(X=list(X=X, model="FIXED", saveEffects=TRUE), 
            G=list(X=W, model="BRR", saveEffects=TRUE),
            E=list(X=E, model="BRR", saveEffects=TRUE)
)

### Ensure ETA$X$X is numeric ###
if (!is.numeric(ETA$X$X)) {
  ETA$X$X <- as.matrix(ETA$X$X)
  ETA$X$X <- apply(ETA$X$X, 2, as.numeric)
}

print("model ETA created")

### Collect results ###

### Model parameters ###  

### intercept ###
zz0 <- read.table(paste(scratch, '/PP_plink_run_G_E_mu.dat', sep=''), header=F); colnames(zz0) <- "int"

print("zz0 created")

### variance of ridge regression betas for the additive genetic effect ###
zz1 <- read.table(paste(scratch, '/PP_plink_run_G_E_ETA_G_varB.dat', sep=''), header=F); colnames(zz1) <- "G"

print("zz1 created")

### variance of ridge regression betas for the environmental effect ###
zz2 <- read.table(paste(scratch, '/PP_plink_run_G_E_ETA_E_varB.dat', sep=''), header=F); colnames(zz2) <- "E"

### variance of the residual deviations ###
zz9 <- read.table(paste(scratch, '/PP_plink_run_G_E_varE.dat', sep=''), header=F); colnames(zz9) <- "res"

print("zz9 created")

### dataframe with model parameters ###
VCEm_PP_plink <- data.frame(zz0, zz1, zz2, zz9)

print("VCEm_PP_plink created now saving")

write.csv(VCEm_PP_plink, file="/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/VCEm_PP_plink_run_G_E.csv", row.names=TRUE)
print("VCEm_PP_plink written")

### Sampled regression effects ###  

### intercept and fixed effects ###
B1 <- read.table(paste(scratch, '/PP_plink_run_G_E_ETA_X_b.dat', sep=''), header=T)

print("B1 done now B2")

### additive genetic random effects ###
B2 <- readBinMat(paste(scratch, '/PP_plink_run_G_E_ETA_G_b.bin', sep=''))

print("B2 done")

### environmental random effects ###
B3 <- readBinMat(paste(scratch, '/PP_plink_run_G_E_ETA_E_b.bin', sep=''))

### dataframe with variance partition ###
varabs <- matrix(NA, nrow_varabs, 3); colnames(varabs) <- c("V_X", "V_G", "V_E")

print("filling up cols of varab")

### fill up column for fixed covariates ###
varabs[, 1] <- matrixStats::colVars(ETA$X$X%*%t(B1))[-c(1:(burnin/thin))]   
print("varab col one filled up now col 2")

### fill up column for random additive genetic covariates ###
varabs[, 2] <- matrixStats::colVars(tcrossprod(ETA$G$X,B2))

### fill up column for random environmental covariates ###
varabs[, 3] <- matrixStats::colVars(tcrossprod(ETA$E$X,B3))

print("varab cols done now saving varab")

write.csv(varabs, file="/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_PP_plink_run_G_E.csv", row.names=TRUE)

print("varab saved")

### Note -c(1:(burnin/thin)). This is because burn-in samples are stored in B1, but are not in B2 and B3.