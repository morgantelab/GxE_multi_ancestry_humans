rm(list=ls()); gc()

library(matrixcalc)
library(genio)
library(data.table)
library(Matrix)

# Load precomputed G and E matrices
#load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/")  # G
plink_grm_path <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/grm_by_plink"
grm_data <- read_grm(plink_grm_path)
load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/Emat_20241106.RData")  # E

#G <- grm_matrix_pcrelate_5pcs
G <- grm_data
G <- G$kinship
E <- Emat

# Ensure G and E are numeric matrices
if (!is.matrix(G)) G <- as.matrix(G)
if (!is.numeric(G)) G <- as.numeric(G)

E <- tcrossprod(Emat)
E <- E/mean(diag(E))

if (!is.matrix(E)) E <- as.matrix(E)
if (!is.numeric(E)) E <- as.numeric(E)

# Check dimensions match
if (nrow(G) != nrow(E) || ncol(G) != ncol(E)) {
  stop("Error: Dimensions of G and E do not match.")
}

# Calculate the GxE interaction matrix
GE <- hadamard.prod(G, E)  # Element-wise multiplication of G and E

# Save GxE results
save(GE, file = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/GE_hadamard_prod_plink.RData")
