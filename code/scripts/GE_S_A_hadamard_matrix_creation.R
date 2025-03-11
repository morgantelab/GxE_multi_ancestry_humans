rm(list=ls()); gc()
set.seed(1123)

library(matrixcalc)
library(genio)
library(data.table)
library(Matrix)

#load G of pcrelate
load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/grm_matrix_pcrelate_5pcs.RData")
G_pcrelate <- grm_matrix_pcrelate_5pcs

if (!is.matrix(G_pcrelate)) G_pcrelate <- as.matrix(G_pcrelate)
if (!is.numeric(G_pcrelate)) G <- as.numeric(G_pcrelate)

# load Emat with sex and age
E <- readRDS("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/Emat_with_sex_age.RDS")
E <- tcrossprod(E)
E <- E/mean(diag(E))

if (!is.matrix(E)) E <- as.matrix(E)
if (!is.numeric(E)) E <- as.numeric(E)

# check dimensions
if (nrow(G_pcrelate) != nrow(E) || ncol(G_pcrelate) != ncol(E)) {
  stop("Error: Dimensions of G and E do not match.")
}

# create GE hadamard product matrix
GE_pcrelate <- hadamard.prod(G_pcrelate, E)  # Element-wise multiplication of G and E

# eigen of GE_pcrelate
GE_pcrelate_eigen <- eigen(GE_pcrelate)

# Assuming row names are already set as individual IDs from 'G'
rownames(GE_pcrelate_eigen$vectors) <- rownames(GE_pcrelate)

# save GE matrix
save(GE_pcrelate, file = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/GE_S_A_hadamard_prod_pcrelate.RDS")

#save
saveRDS(GE_pcrelate_eigen, file = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/GE_S_A_eigen_pcrelate.RDS")

