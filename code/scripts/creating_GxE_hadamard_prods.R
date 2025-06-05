rm(list=ls()); gc()

library(matrixcalc)
library(genio)
library(data.table)
library(Matrix)

# Load precomputed G and E matrices
#load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/")  # G
plink_grm_path <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/grm_by_plink"
plink_grm_data <- read_grm(plink_grm_path)
load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/Emat_20250106.RData")  # E

# Load precomputed G and E matrices
#load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/")  # G
load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/grm_matrix_pcrelate_5pcs.RData")
G_pcrelate <- grm_matrix_pcrelate_5pcs
G_plink <- plink_grm_data
G_plink <- G_plink$kinship
E <- Emat

# Ensure G and E are numeric matrices
if (!is.matrix(G_plink)) G_plink <- as.matrix(G_plink)
if (!is.numeric(G_plink)) G <- as.numeric(G_plink)
if (!is.matrix(G_pcrelate)) G_pcrelate <- as.matrix(G_pcrelate)
if (!is.numeric(G_pcrelate)) G <- as.numeric(G_pcrelate)


E <- tcrossprod(Emat)
E <- E/mean(diag(E))

if (!is.matrix(E)) E <- as.matrix(E)
if (!is.numeric(E)) E <- as.numeric(E)

# Check dimensions match
if (nrow(G_plink) != nrow(E) || ncol(G_plink) != ncol(E)) {
  stop("Error: Dimensions of G and E do not match.")
}

if (nrow(G_pcrelate) != nrow(E) || ncol(G_pcrelate) != ncol(E)) {
  stop("Error: Dimensions of G and E do not match.")
}

# Calculate the GxE interaction matrix
GE_plink <- hadamard.prod(G_plink, E)  # Element-wise multiplication of G and E
GE_pcrelate <- hadamard.prod(G_pcrelate, E)  # Element-wise multiplication of G and E

# Save GxE results
save(GE_plink, file = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/GE_hadamard_prod_plink.RData")
save(GE_pcrelate, file = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/GE_hadamard_prod_pcrelate.RData")
