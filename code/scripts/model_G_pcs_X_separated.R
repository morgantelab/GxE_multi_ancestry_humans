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
  make_option(c("-p", "--pcs"), type = "character", default = NULL, 
              help = "Path to the scaled PCs RDS file", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL, 
              help = "Output directory for saving results", metavar = "character"),
  make_option(c("-s", "--scratch"), type = "character", default = "/scratch3/kgoda/ukbiobank_files/tmp/snakemake_runs", 
              help = "Temporary directory for storing model files", metavar = "character"),
  make_option(c("-r", "--result"), type = "character", default = NULL, 
              help = "taking in output file name to get type", metavar = "character")
  
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check required arguments
if (is.null(opt$data) | is.null(opt$eigen) | is.null(opt$pcs) | is.null(opt$output)) {
  print_help(opt_parser)
  stop("All file paths and phenotype type must be provided.", call. = FALSE)
}

# Set the working directory dynamically
setwd(opt$dir)

# Extract `type` and `grm_source` from output file path
output_file <- basename(opt$result)  # Get the file name from the path
parsed_info <- strsplit(output_file, "_")[[1]]  # Split the file name by underscore

# Assuming the file name follows the pattern: "varabs_<type>_<grm_source>_G_pcs.csv"
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
X <- matched_dataset[, c("AOPs", "AOPss", "Sex_SI")]
X$AOPs <- as.vector(X$AOPs)
X$AOPss <- as.vector(X$AOPss)
rownames(X) <- rownames(W)

# Load scaled PCs and match to individual IDs
pcs_scaled <- readRDS(opt$pcs)
matched_pcs <- pcs_scaled[match(individual_ids, pcs_scaled$ID), 2:11]
P <- matched_pcs

print("Combined X matrix created")

# Model setup
iter <- 90000
burnin <- 40000
thin <- 50
verb <- T
nrow_varabs <- (iter-burnin)/thin

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

model <- BGLR(y=y, ETA=ETA, nIter=iter, burnIn=burnin, thin=thin, verbose=verb, saveAt=paste(opt$scratch, '/', type, '_', grm_source, '_run_G_pcs_X_separated_', sep=''))

# Collect results
#zz0 <- read.table(paste(opt$scratch, '/', type, '_', grm_source, '_run_G_pcs_mu.dat', sep=''), header=F)
#colnames(zz0) <- "int"
#zz1 <- read.table(paste(opt$scratch, '/', type, '_', grm_source, '_run_G_pcs_ETA_G_varB.dat', sep=''), header=F)
#colnames(zz1) <- "G"
#zz9 <- read.table(paste(opt$scratch, '/', type, '_', grm_source, '_run_G_pcs_varE.dat', sep=''), header=F)
#colnames(zz9) <- "res"

#VCEm <- data.frame(zz0, zz1, zz9)
#write.csv(VCEm, file=file.path(opt$output, paste0("VCEm_", type, "_", grm_source, "_G_pcs.csv")), row.names=TRUE)
#print("VCEm results written")

# Save variance partition results
# intercept and fixed effects
#B1 <- read.table(paste(opt$scratch, '/', type, '_', grm_source, '_run_G_pcs_ETA_X_B.dat', sep=''), header=T)

# additive genetic random effects
#B2 <- readBinMat(paste(opt$scratch, '/', type, '_', grm_source, '_run_G_pcs_ETA_G_b.bin', sep=''))

#varabs <- matrix(NA, nrow_varabs, 2); colnames(varabs) <- c("V_X", "V_G")
#varabs[, 1] <- matrixStats::colVars(ETA$X$X %*% t(B1))[-c(1:(burnin/thin))]
#varabs[, 2] <- matrixStats::colVars(tcrossprod(ETA$G$X, B2))
#write.csv(varabs, file=file.path(opt$output, paste0("varabs_", type, "_", grm_source, "_G_pcs.csv")), row.names=TRUE)

