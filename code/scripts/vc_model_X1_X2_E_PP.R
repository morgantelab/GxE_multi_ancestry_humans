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
load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/scaled_dataset_20250106.Rdata")
print("dataset loaded")

### Extract the phenotype vectors ###
y <- dataset$PP0s
rownames(y) <- dataset$ID

print("y_PP_pcrelate created")

### X is the incidence matrix for the 'fixed' covariates (no penalisation, no shrinkage). here age and sex ###
### Extract the covariate matrix ###
X <- dataset[, c("AOPs", "AOPss", "Sex_SI")]
X$AOPs <- as.vector(X$AOPs)
X$AOPss <- as.vector(X$AOPss)
rownames(X) <- dataset$ID

print("initial X created")

# Load scaled PCs and match to individual IDs
pcs_scaled <- readRDS("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/scaled_pcs_plink.rds")
matched_pcs <- pcs_scaled[match(dataset$ID, pcs_scaled$ID), 2:11]
P <- matched_pcs

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


iter <- 90000
burnin <- 40000
thin <- 50
verb <- T
nrow_varabs <- (iter-burnin)/thin
scratch <- "/scratch3/kgoda/ukbiobank_files/tmp/snakemake_runs" # Use this!

### Define model ###

ETA <- list(X1=list(X=X, model="FIXED", saveEffects=TRUE),
            X2=list(X=P, model="FIXED", saveEffects=TRUE),
            E=list(X=E, model="BRR", saveEffects=TRUE)
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

model <- BGLR(y=y, ETA=ETA, nIter=iter, burnIn=burnin, thin=thin, verbose=verb, saveAt=paste(scratch, '/PP_X1_X2_E_run_', sep=''))

# # ### Collect results ###
# #
# # ### Model parameters ###
# #
# # ### intercept ###
# # zz0 <- read.table(paste(scratch, '/PP_pcrelate_run_G_E_mu.dat', sep=''), header=F); colnames(zz0) <- "int"
# #
# # print("zz0 created")
# #
# # ### variance of ridge regression betas for the additive genetic effect ###
# # zz1 <- read.table(paste(scratch, '/PP_pcrelate_run_G_E_ETA_G_varB.dat', sep=''), header=F); colnames(zz1) <- "G"
# #
# # print("zz1 created")
# #
# # ### variance of ridge regression betas for the environmental effect ###
# # zz2 <- read.table(paste(scratch, '/PP_pcrelate_run_G_E_ETA_E_varB.dat', sep=''), header=F); colnames(zz2) <- "E"
# #
# # ### variance of the residual deviations ###
# # zz9 <- read.table(paste(scratch, '/PP_pcrelate_run_G_E_varE.dat', sep=''), header=F); colnames(zz9) <- "res"
# #
# # print("zz9 created")
# #
# # ### dataframe with model parameters ###
# # VCEm_PP_pcrelate <- data.frame(zz0, zz1, zz2, zz9)
# #
# # print("VCEm_PP_pcrelate created now saving")
# #
# # write.csv(VCEm_PP_pcrelate, file="/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/VCEm_PP_pcrelate_run_G_E.csv", row.names=TRUE)
# # print("VCEm_PP_pcrelate written")
# #
# # ### Sampled regression effects ###
# #
# # ### intercept and fixed effects ###
# ### Sampled regression effects ###
#
# ### intercept and fixed effects ###
# B1 <- read.table(paste(scratch, '/PP_X1_run_ETA_X_b.dat', sep=''), header=T)
#
# ### dataframe with variance partition ###
# varabs <- matrix(NA, nrow_varabs, 1); colnames(varabs) <- c("V_X")
#
# print("filling up cols of varab")
#
# ### fill up column for fixed covariates ###
# varabs[, 1] <- matrixStats::colVars(ETA$X$X%*%t(B1))[-c(1:(burnin/thin))]
#
# print("varab cols done now saving varab")
#
# write.csv(varabs, file="/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_PP_X.csv", row.names=TRUE)
#
# print("varab saved")
#
#
# # ### Note -c(1:(burnin/thin)). This is because burn-in samples are stored in B1, but are not in B2 and B3.
