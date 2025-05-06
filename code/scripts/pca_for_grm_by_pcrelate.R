set.seed(1123)

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

print("loading grm")
# Load GRM matrix from the provided input file
readRDS(opt$input)
print("grm loaded")

print("computing eigen")
# Compute eigenvalues and eigenvectors
eigen_results <- eigen(grm_matrix_pcrelate_5pcs)
print("eigen calculated")

print("fixing rownames")
# Assuming row names are already set as individual IDs from the GRM matrix
rownames(eigen_results$vectors) <- rownames(grm_matrix_pcrelate_5pcs)

print("saving eigen results")
# Save the eigenvalues and eigenvectors to the output RDS file
saveRDS(eigen_results, file = opt$output)
print("eigen saved")

