# Load necessary libraries
library(genio)
library(data.table)
library(Matrix)
library(optparse)

# Define command-line options
option_list <- list(
  make_option(c("-d", "--dir"), type = "character", default = NULL, 
              help = "path to the working directory", metavar = "character"),
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "path to the GRM matrix RData file", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "path to save the eigenvalue and eigenvector results", metavar = "character")
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if the required arguments are provided
if (is.null(opt$input) | is.null(opt$output)) {
  print_help(opt_parser)
  stop("Both input and output file paths must be provided.", call. = FALSE)
}

# Set the working directory dynamically
setwd(opt$dir)

# Read the GRM using genio
grm_data <- read_grm(opt$input)

# Perform Eigen decomposition (PCA) on the kinship matrix
eigen_results <- eigen(grm_data$kinship)

# Set individual IDs as rownames for the eigenvectors
rownames(eigen_results$vectors) <- grm_data$fam$id

# Save the eigenvalues and eigenvectors to an RDS file
saveRDS(eigen_results, file = opt$output)

