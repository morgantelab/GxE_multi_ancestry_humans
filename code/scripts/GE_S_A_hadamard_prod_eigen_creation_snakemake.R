#### Sex and Age added to Emat, GE creation, Eigen of GE. ####

rm(list=ls()); gc()

# Set seed for reproducibility
set.seed(1123)

# Load necessary libraries
library(optparse)
library(matrixcalc)
library(genio)
library(data.table)
library(Matrix)

# Define command-line options
option_list <- list(
  make_option(c("-d", "--outdir"), type="character", default="output/",
              help="Directory where output files will be saved", metavar="character"),
  make_option(c("-g", "--grm"), type="character", default=NULL, 
              help="Path to GRM matrix file", metavar="character"),
  make_option(c("-e", "--emat"), type="character", default=NULL, 
              help="Path to Emat file", metavar="character")
)

# Parse options
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$grm) || is.null(opt$emat)) {
  print_help(opt_parser)
  stop("Missing required arguments.")
}

# Define output file paths
GE_path <- file.path(opt$outdir, "GE_S_A_hadamard_prod_pcrelate.RDS")
GE_eigen_path <- file.path(opt$outdir, "GE_S_A_eigen_pcrelate.RDS")

# Load G matrix
load(opt$grm)
G_pcrelate <- grm_matrix_pcrelate_5pcs

if (!is.matrix(G_pcrelate)) G_pcrelate <- as.matrix(G_pcrelate)
if (!is.numeric(G_pcrelate)) G_pcrelate <- as.numeric(G_pcrelate)

# Load Emat matrix
E <- readRDS(opt$emat)
E <- tcrossprod(E)
E <- E / mean(diag(E))

if (!is.matrix(E)) E <- as.matrix(E)
if (!is.numeric(E)) E <- as.numeric(E)

# Check dimensions
if (nrow(G_pcrelate) != nrow(E) || ncol(G_pcrelate) != ncol(E)) {
  stop("Error: Dimensions of G and E do not match.")
}

# Compute Hadamard product
GE_pcrelate <- hadamard.prod(G_pcrelate, E)

# Compute eigen decomposition
GE_pcrelate_eigen <- eigen(GE_pcrelate)

# Set row names
rownames(GE_pcrelate_eigen$vectors) <- rownames(G_pcrelate)

# Save outputs
saveRDS(GE_pcrelate, file=GE_path)
saveRDS(GE_pcrelate_eigen, file=GE_eigen_path)

