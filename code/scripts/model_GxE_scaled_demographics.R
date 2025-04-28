set.seed(1123)

# Load required libraries
library(data.table)
library(Matrix)
library(BGLR)
library(readr)
library(genio)
library(dplyr)
library(optparse)

# Command-line options
option_list <- list(
  make_option(c("-d", "--dir"), type = "character", default = NULL,
              help = "path to the working directory", metavar = "character"),
  make_option(c("-t", "--data"), type = "character", default = NULL,
              help = "Path to the scaled dataset RData file", metavar = "character"),
  make_option(c("-e", "--eigen"), type = "character", default = NULL,
              help = "Path to the eigen results RDS file", metavar = "character"),
  make_option(c("-v", "--envvars"), type = "character", default = NULL,
              help = "envvar matrix", metavar = "character"),
  make_option(c("-p", "--pcs"), type = "character", default = NULL,
              help = "Path to the scaled PCs RDS file", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output directory for saving results", metavar = "character"),
  make_option(c("-q", "--GE"), type = "character", default = NULL,
              help = "Path to the GE hadamad pdt", metavar = "character"),
  make_option(c("-s", "--scratch"), type = "character", default = "/scratch3/kgoda/ukbiobank_files/tmp/snakemake_runs",
              help = "Temporary directory for storing model files", metavar = "character"),
  make_option(c("-r", "--result"), type = "character", default = NULL,
              help = "taking in output file name to get type", metavar = "character")

)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check required arguments
if (is.null(opt$data) | is.null(opt$eigen) | is.null(opt$scratch) | is.null(opt$output)) {
  print_help(opt_parser)
  stop("All file paths and phenotype type must be provided.", call. = FALSE)
}

# Set the working directory dynamically
setwd(opt$dir)

# Extract `type` and `grm_source` from output file path
output_file <- basename(opt$result)  # Get the file name from the path
parsed_info <- strsplit(output_file, "_")[[1]]  # Split the file name by underscore

# Assuming the file name follows the pattern: "varabs_<type>_<grm_source>_GE.csv"
type <- parsed_info[2]  # Extract the type (e.g., sp, dp, pp)
grm_source <- parsed_info[3]  # Extract the GRM source (e.g., plink, pcrelate)

# Load dataset
load(opt$data)

# Load eigen results
eigen_results <- readRDS(opt$eigen)
eigenvectors <- eigen_results$vectors
eigenvalues <- eigen_results$values

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

# Prepare phenotype vector based on type
y <- matched_dataset[[paste0(type, "0s")]]

# Prepare covariate matrix
X <- matched_dataset[, c("AOPs", "AOPsss", "Sex_SIs")]
X$AOPs <- as.vector(X$AOPs)
X$AOPss <- as.vector(X$AOPsss)
rownames(X) <- rownames(W)

# Load scaled PCs and match to individual IDs
pcs_scaled <- readRDS(opt$pcs)
matched_pcs <- pcs_scaled[match(individual_ids, pcs_scaled$ID), 2:11]
P <- matched_pcs

# Load eigen of E
E_eigen <- readRDS(opt$envvars)
E_eigenvectors <- E_eigen$vectors
E_eigenvalues <- E_eigen$values

# Filter and scale eigenvectors by positive eigenvalues
E_positive_indices <- which(E_eigenvalues > 0)
E_filtered_eigenvectors <- E_eigenvectors[, E_positive_indices]
E_filtered_eigenvalues <- E_eigenvalues[E_positive_indices]
for(i in 1:ncol(E_filtered_eigenvectors)) {
  E_filtered_eigenvectors[, i] <- E_filtered_eigenvectors[, i] * sqrt(E_filtered_eigenvalues[i])
}
E <- E_filtered_eigenvectors

# Load eigen of GE
GE_eigen <- readRDS(opt$GE)
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

ETA <- list(X1=list(X=X, model="FIXED", saveEffects=TRUE),
            X2=list(X=P, model="FIXED", saveEffects=TRUE),
            G=list(X=W, model="BRR", saveEffects=TRUE),
            E=list(X=E, model="BRR", saveEffects=TRUE),
            GxE=list(X=GE, model="BRR", saveEffects=TRUE)
)

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

model <- BGLR(y=y, ETA=ETA, nIter=iter, burnIn=burnin, thin=thin, verbose=verb, saveAt=paste(opt$scratch, '/', type, '_', grm_source, '_run_GxE_pcrelate_pcs_plink_scaled_demographics_', sep=''))

###Collecting results###
# Load results from BGLR output
#zz0 <- read.table(paste(opt$scratch, '/', type, '_', grm_source, '_run_GxE_pcrelate_pcs_plink_scaled_demographics_mu.dat', sep=''), header=F)
#colnames(zz0) <- "int"

#zz1 <- read.table(paste(opt$scratch, '/', type, '_', grm_source, '_run_GxE_pcrelate_pcs_plink_scaled_demographics_ETA_G_varB.dat', sep=''), header=F)
#colnames(zz1) <- "G"

#zz2 <- read.table(paste(opt$scratch, '/', type, '_', grm_source, '_run_GxE_pcrelate_pcs_plink_scaled_demographics_ETA_E_varB.dat', sep=''), header=F)
#colnames(zz2) <- "E"

#zz9 <- read.table(paste(opt$scratch, '/', type, '_', grm_source, '_run_GxE_pcrelate_pcs_plink_scaled_demographics_varE.dat', sep=''), header=F)
#colnames(zz9) <- "res"

# Combine results into a dataframe and save
#VCEm <- data.frame(zz0, zz1, zz2, zz9)
#write.csv(VCEm, file=file.path(opt$output, paste0("VCEm_", type, "_", grm_source, "_run_GxE_pcrelate_pcs_plink_scaled_demographics.csv")), row.names=TRUE)

# Sampled regression effects
B1 <- read.table(paste(opt$scratch, '/', type, '_', grm_source, '_run_GxE_pcrelate_pcs_plink_scaled_demographics_ETA_X1_b.dat', sep=''), header=TRUE)
B2 <- read.table(paste(opt$scratch, '/', type, '_', grm_source, '_run_GxE_pcrelate_pcs_plink_scaled_demographics_ETA_X2_b.dat', sep=''), header=TRUE)
B3 <- readBinMat(paste(opt$scratch, '/', type, '_', grm_source, '_run_GxE_pcrelate_pcs_plink_scaled_demographics_ETA_G_b.bin', sep=''))
B4 <- readBinMat(paste(opt$scratch, '/', type, '_', grm_source, '_run_GxE_pcrelate_pcs_plink_scaled_demographics_ETA_E_b.bin', sep=''))
B5 <- readBinMat(paste(opt$scratch, '/', type, '_', grm_source, '_run_GxE_pcrelate_pcs_plink_scaled_demographics_ETA_GxE_b.bin', sep=''))

# Calculate variance components
varabs <- matrix(NA, nrow_varabs, 5); colnames(varabs) <- c("V_X1", "V_X2", "V_G", "V_E", "V_GxE")

# Fill variance components
varabs[, 1] <- matrixStats::colVars(ETA$X1$X %*% t(B1))[-c(1:(burnin/thin))]
varabs[, 2] <- matrixStats::colVars(ETA$X2$X %*% t(B2))[-c(1:(burnin/thin))]
varabs[, 3] <- matrixStats::colVars(tcrossprod(ETA$G$X, B3))
varabs[, 4] <- matrixStats::colVars(tcrossprod(ETA$E$X, B4))
varabs[, 5] <- matrixStats::colVars(tcrossprod(ETA$GxE$X, B5))

# Save variance components
write.csv(varabs, file=file.path(opt$output, paste0("varabs_", type, "_", grm_source, "_GxE_pcrelate_pcs_plink_scaled_demographics.csv")), row.names=TRUE)
print("Variance partition results saved")
