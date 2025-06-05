set.seed(1123)

# Clear environment
rm(list=ls()); gc()

# Load libraries
library(matrixcalc)
library(genio)
library(data.table)
library(Matrix)
library(optparse)

# Define options
option_list <- list(
  make_option(c("--plink_grm"), type = "character", help = "Path to PLINK GRM prefix"),
  make_option(c("--pcrelate_grm"), type = "character", help = "Path to PCRelate GRM RDS file"),
  make_option(c("--emat"), type = "character", help = "Path to scaled E matrix RData file"),
  make_option(c("--out_GxE_plink"), type = "character", help = "Output path for GE_plink RData"),
  make_option(c("--out_GxE_pcrelate"), type = "character", help = "Output path for GE_pcrelate RData"),
  make_option(c("--out_eigen_plink"), type = "character", help = "Output path for eigen GE_plink RDS"),
  make_option(c("--out_eigen_pcrelate"), type = "character", help = "Output path for eigen GE_pcrelate RDS")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Load GRMs and E matrix
plink_grm_data <- read_grm(opt$plink_grm)
G_plink <- plink_grm_data$kinship
G_pcrelate <- readRDS(opt$pcrelate_grm)
load(opt$emat)

# Compute E from Emat
E <- tcrossprod(Emat)
E <- E / mean(diag(E))

# Ensure numeric matrices
G_plink <- as.matrix(G_plink)
G_pcrelate <- as.matrix(G_pcrelate)
E <- as.matrix(E)

# Dimension checks
if (!all(dim(G_plink) == dim(E))) stop("Mismatch in G_plink and E dimensions")
if (!all(dim(G_pcrelate) == dim(E))) stop("Mismatch in G_pcrelate and E dimensions")

# Hadamard product: GxE
GE_plink <- hadamard.prod(G_plink, E)
GE_pcrelate <- hadamard.prod(G_pcrelate, E)

# Save GE matrices
save(GE_plink, file = opt$out_GxE_plink)
save(GE_pcrelate, file = opt$out_GxE_pcrelate)

# Eigen decompositions
GE_plink_eigen <- eigen(GE_plink)
rownames(GE_plink_eigen$vectors) <- rownames(GE_plink)
saveRDS(GE_plink_eigen, file = opt$out_eigen_plink)

GE_pcrelate_eigen <- eigen(GE_pcrelate)
rownames(GE_pcrelate_eigen$vectors) <- rownames(GE_pcrelate)
saveRDS(GE_pcrelate_eigen, file = opt$out_eigen_pcrelate)
