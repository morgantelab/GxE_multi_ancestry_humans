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
              help = "Path to the scaled dataset RDS file", metavar = "character"),
  make_option(c("-p", "--pcs"), type = "character", default = NULL,
              help = "Output directory for saving results", metavar = "character"),
  make_option(c("-e", "--eigen"), type = "character", default = NULL,
              help = "Path to the eigen results RDS file", metavar = "character"),
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
if (is.null(opt$data) | is.null(opt$output)) {
  print_help(opt_parser)
  stop("All file paths and phenotype type must be provided.", call. = FALSE)
}

# Set the working directory dynamically
setwd(opt$dir)

# Extract `type` and `grm_source` from output file path
output_file <- basename(opt$result)  # Get the file name from the path
parsed_info <- strsplit(output_file, "_")[[1]]  # Split the file name by underscore

# Assuming the file name follows the pattern: "VCEm_<type>_<grm_source>_G.csv"
type <- parsed_info[2]  # Extract the type (e.g., sp, dp, pp)
grm <- parsed_info[3]  # Extract the GRM source (e.g., plink, pcrelate)

# Load dataset
dataset<- readRDS(opt$data)

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
rownames(y) <- rownames(W)

# Prepare covariate matrix
X <- matched_dataset[, c("AOPs", "AOPsss", "Sex_SIs")]
X$AOPs <- as.vector(X$AOPs)
X$AOPsss <- as.vector(X$AOPsss)
X$Sex_SIs <- as.vector(X$Sex_SIs)
rownames(X) <- rownames(W)

# Load scaled PCs and match to individual IDs
pcs_scaled <- readRDS(opt$pcs)
P <- pcs_scaled[match(individual_ids, pcs_scaled$ID), 2:11]

# Model setup
iter <- 90000
burnin <- 40000
thin <- 50
verb <- T
nrow_varabs <- (iter-burnin)/thin

ETA <- list(X1=list(X=X, model="FIXED", saveEffects=TRUE),
            X2=list(X=P, model="FIXED", saveEffects=TRUE),
            G=list(X=W, model="BRR", saveEffects=TRUE)
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
print("Model ETA created")

# Run BGLR model

model <- BGLR(y=y, ETA=ETA, nIter=iter, burnIn=burnin, thin=thin, verbose=verb, saveAt=paste(opt$scratch, '/', type, '_', grm, '_X1_X2_G_', sep=''))

## Collecting results ##

### Sampled regression effects ###

# Construct the correct file path dynamically
model_output_prefix <- paste0(opt$scratch, "/", type, '_', grm, "_X1_X2_G_")

### additive genetic random effects ###
B1 <- read.table(paste0(model_output_prefix, "ETA_X1_b.dat"), header = TRUE)
B2 <- read.table(paste0(model_output_prefix, "ETA_X2_b.dat"), header = TRUE)
B3 <- readBinMat(paste0(model_output_prefix, "ETA_G_b.bin"))

### dataframe with variance partition ###
varabs <- matrix(NA, nrow_varabs, 3); colnames(varabs) <- c("V_X1", "V_X2", "V_G")

print("filling up cols of varab")

# Fill variance components
varabs[, 1] <- matrixStats::colVars(ETA$X1$X %*% t(B1))[-c(1:(burnin/thin))]
varabs[, 2] <- matrixStats::colVars(ETA$X2$X %*% t(B2))[-c(1:(burnin/thin))]
varabs[, 3] <- matrixStats::colVars(tcrossprod(ETA$G$X, B3))

print("varab cols done now saving varab")

# Dynamic output filename
output_file <- file.path(opt$output, paste0("varabs_", type, '_', grm, "_X1_X2_G", ".csv"))

# Save variance components
write.csv(varabs, file = output_file, row.names = TRUE)

print("varab saved")
