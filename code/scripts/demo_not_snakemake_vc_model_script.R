### For PLINK and PCRELATE GRM ###
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
eigen_path_pcrelate <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/pca_for_pcrelate.rds"
eigen_results_pcrelate <- readRDS(eigen_path_pcrelate)
eigenvectors <- eigen_results_pcrelate$vectors
eigenvalues <- eigen_results_pcrelate$values
print("eigen results handled")

### For GBLUP, W is a matrix of eigenvectors each multiplied by the respective eigenvalues. This implies that the first eigenvectors will have larger variance. Also, please remove eigenvectors with negative eigenvalues.
### Fabio recommended using this to create W ###
### W=eigen_results$vectors
### for(i in 1:ncol(W)){  W[,i]=W[,i]*sqrt(eigen_results$values[i]) }
### W=W[,eigen_results$values>1e-5]

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
y_PP_pcrelate <- matched_dataset$PP0s

print("y_PP_pcrelate created")

### X is the incidence matrix for the 'fixed' covariates (no penalisation, no shrinkage). here age and sex ###
### Extract the covariate matrix ###
X <- matched_dataset[, c("AOPs", "AOPss", "Sex_SI")]
X$AOPs <- as.vector(X$AOPs)
X$AOPss <- as.vector(X$AOPss)
rownames(X) <- rownames(W)

print("initial X created")

#load scaled PCs
#pcs_scaled <- readRDS("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/scaled_pcs_plink.rds")
#print("scaled pcs loaded")

# match the PCs to the individual_ids
#matched_pcs <- pcs_scaled[match(individual_ids, pcs_scaled$ID), ]

# Extract only the PC columns
#pcs <- matched_pcs[, 2:11]  # Assuming the first column is ID and the next 10 columns are PC1-PC10

# Combine the covariates and PCs
#X <- cbind(X, pcs)

# Print the structure of X to verify
#print("Combined X matrix created")
#str(X)

load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/Emat_20241106.RData")
filtered_Emat <- Emat[match(individual_ids, rownames(Emat)), ]
rownames(filtered_Emat) <- individual_ids
E <- tcrossprod(filtered_Emat) / ncol(filtered_Emat)

### small checks for the matrices and vectors created ###
#dataset[dataset$ID == 1000086, "PP0a"]
#dataset[dataset$ID == 1000086, "DP0a"]
#dataset[dataset$ID == 1000086, "DP0a"]
#dataset[dataset$ID == 1000086, "AOP"]
#dataset[dataset$ID == 1000086, "Sex_SI"]

### Model ###

iter <- 90000
burnin <- 40000
thin <- 50
verb <- T
nrow_varabs <- (iter-burnin)/thin
#nrow_varabs <- 1000
scratch <- "/scratch3/kgoda/ukbiobank_files/tmp/snakemake_runs" # Use this!

### Define model ###

ETA <- list(X=list(X=X, model="FIXED", saveEffects=TRUE), 
            G=list(X=W, model="BRR", saveEffects=TRUE)
)

### Ensure ETA$X$X is numeric ###
if (!is.numeric(ETA$X$X)) {
  ETA$X$X <- as.matrix(ETA$X$X)
  ETA$X$X <- apply(ETA$X$X, 2, as.numeric)
}

print("model ETA created")

### Run model here ###
#print("model for PP_pcrelate running")

#model <- BGLR(y=y_PP_pcrelate, ETA=ETA, nIter=iter, burnIn=burnin, thin=thin, verbose=verb, saveAt=paste(scratch, '/PP_pcrelate_G_run_', sep=''))

#print("model for PP_pcrelate done")

### Collect results ###

### Model parameters ###  

### intercept ###
zz0 <- read.table(paste(scratch, '/PP_pcrelate_G_run_mu.dat', sep=''), header=F); colnames(zz0) <- "int"

print("zz0 created")

### variance of ridge regression betas for the additive genetic effect ###
zz1 <- read.table(paste(scratch, '/PP_pcrelate_G_run_ETA_G_varB.dat', sep=''), header=F); colnames(zz1) <- "G"

print("zz1 created")

### variance of ridge regression betas for the environmental effect ###
#zz2 <- read.table(paste(scratch, '/PP_pcrelate_run_ETA_E_varB.dat', sep=''), header=F); colnames(zz2) <- "E"

### variance of the residual deviations ###
zz9 <- read.table(paste(scratch, '/PP_pcrelate_G_run_varE.dat', sep=''), header=F); colnames(zz9) <- "res"

print("zz9 created")

### dataframe with model parameters ###
VCEm_PP_pcrelate <- data.frame(zz0, zz1, zz9)
#VCEm_PP_pcrelate <- data.frame(zz0, zz1, zz9)

print("VCEm_PP_pcrelate created now saving")

write.csv(VCEm_PP_pcrelate, file="/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/VCEm_PP_pcrelate_G.csv", row.names=TRUE)
#VCEm_PP <- read.csv("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/apr_run_with_chinese_20240402/5_variance_components/VCEm_PP_pcrelate")
print("VCEm_PP_pcrelate written")

### Sampled regression effects ###  

### intercept and fixed effects ###
B1 <- read.table(paste(scratch, '/PP_pcrelate_G_run_ETA_X_b.dat', sep=''), header=T)

print("B1 done now B2")

### additive genetic random effects ###
B2 <- readBinMat(paste(scratch, '/PP_pcrelate_G_run_ETA_G_b.bin', sep=''))

print("B2 done")

### environmental random effects ###
#B3 <- readBinMat(paste(scratch, '/PP_pcrelate_r5_run_this_ETA_E_b.bin', sep=''))

### dataframe with variance partition ###
#varabs <- matrix(NA, nrow_varabs, 3); colnames(varabs) <- c("V_X", "V_G", "V_E")
varabs <- matrix(NA, nrow_varabs, 2); colnames(varabs) <- c("V_X", "V_G")

print("filling up cols of varab")

### fill up column for fixed covariates ###
varabs[, 1] <- matrixStats::colVars(ETA$X$X%*%t(B1))[-c(1:(burnin/thin))]   
print("varab col one filled up now col 2")

### fill up column for random additive genetic covariates ###
varabs[, 2] <- matrixStats::colVars(tcrossprod(ETA$G$X,B2))

# Ensure ETA$E$X is a numeric matrix (doesnt run without this)
#ETA$E$X <- as.matrix(ETA$E$X)
#ETA$E$X <- apply(ETA$E$X, 2, as.numeric)

# Ensure B3 is numeric
#B3 <- as.matrix(B3)
#B3 <- apply(B3, 2, as.numeric)

### fill up column for random environmental covariates ###
#varabs[-c(1:(burnin/thin)), 2] <- apply(ETA$E$X%*%t(B3), 2, var)
#varabs[, 3] <- matrixStats::colVars(tcrossprod(ETA$E$X,B3))

print("varab cols done now saving varab")

#plot(1:1000, varabs[,1])
#plot(1:1000, varabs[,2])

write.csv(varabs, file="/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_PP_pcrelate_G.csv", row.names=TRUE)

print("varab saved")

### fill up column for random environmental covariates ###
#varabs[-c(1:(burnin/thin)), 2] <- apply(ETA$E$X%*%t(B3), 2, var)



### Note -c(1:(burnin/thin)). This is because burn-in samples are stored in B1, but are not in B2 and B3.