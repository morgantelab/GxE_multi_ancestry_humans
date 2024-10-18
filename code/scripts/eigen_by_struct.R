# Load necessary libraries
library(genio)
library(data.table)
library(dplyr)
library(optparse)

# Define command-line options
option_list <- list(
  make_option(c("-d", "--dir"), type = "character", default = NULL, 
              help = "path to the working directory", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL, 
              help = "path to the output file", metavar = "character"),
  make_option(c("-b", "--bfile"), type = "character", default = NULL, 
              help = "path to the Plink binary file prefix", metavar = "character"),
  make_option(c("-s", "--dataset"), type = "character", default = NULL, 
              help = "path to the dataset RData file", metavar = "character"),
  make_option(c("-i", "--idfile"), type = "character", default = NULL, 
              help = "path to the subset of selected ids", metavar = "character")
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if the required arguments are provided
if (is.null(opt$dir) | is.null(opt$output) | is.null(opt$bfile) | is.null(opt$dataset) | is.null(opt$idfile)) {
  print_help(opt_parser)
  stop("the working directory, input files, and output file path must be provided", call. = FALSE)
}

# Load the subset of selected IDs (only first column)
selected_ids <- fread(opt$idfile, header = FALSE, col.names = c("ID1", "ID2"))$ID1

# Define file paths dynamically using command-line arguments
genotype_file_prefix <- file.path(opt$bfile)
ethnicity_file_path <- file.path(opt$dataset)

# Load genotype data using genio
plink_data <- read_plink(genotype_file_prefix)

# Extract individual IDs from the plink data (column names after transposition)
individual_ids <- colnames(plink_data$X)

# Subset the genotype data based on the selected IDs
selected_genotype_indices <- match(selected_ids, individual_ids)
selected_genotype_indices <- selected_genotype_indices[!is.na(selected_genotype_indices)]
if (length(selected_genotype_indices) == 0) {
  stop("No matching individuals found in the genotype data for the provided ID list.")
}
plink_data$X <- plink_data$X[, selected_genotype_indices]  # Subset genotype matrix

# Load ethnicity data from an RData file
load(ethnicity_file_path)  # This will load an object into the environment
new_dataset <- get(ls()[1])  # Adjust based on the actual loaded object

# Ensure that the loaded ethnicity data has the right variable names and types
new_dataset$ID <- as.character(new_dataset$ID)  # Convert ID to character for matching

# Subset the ethnicity data based on the selected IDs
new_dataset <- new_dataset[new_dataset$ID %in% selected_ids, ]
if (nrow(new_dataset) == 0) {
  stop("No matching individuals found in the ethnicity data for the provided ID list.")
}

# Function to impute missing values in the genotype matrix
impute_with_row_means <- function(genotypes) {
  for (i in 1:nrow(genotypes)) {
    missing_indices <- is.na(genotypes[i, ])
    if (any(missing_indices)) {
      genotypes[i, missing_indices] <- mean(genotypes[i, !is.na(genotypes[i, ])], na.rm = TRUE)
    }
  }
  return(genotypes)
}

# Apply the imputation function to the genotype matrix
plink_data$X <- impute_with_row_means(plink_data$X)

# Function to prepare, transpose, and scale genotype matrices by ethnicity
prepare_scale_genotypes <- function(genotypes, target_ethnicity, ethnicity_data) {
    target_ids <- ethnicity_data$ID[ethnicity_data$ethn1_consolidated == target_ethnicity]
    subset_indices <- match(target_ids, colnames(genotypes))
    subset_indices <- subset_indices[!is.na(subset_indices)]

    if (length(subset_indices) == 0) {
        warning(paste("No data found for ethnicity:", target_ethnicity))
        return(NULL)
    }

    # Transpose and scale the imputed genotype matrix
    transposed_genotypes <- t(genotypes[, subset_indices])
    return(scale(transposed_genotypes, center = TRUE, scale = TRUE))
}

# Prepare and scale matrices for each ethnicity, now with imputation
W1 <- prepare_scale_genotypes(plink_data$X, "White", new_dataset)
W2 <- prepare_scale_genotypes(plink_data$X, "Black", new_dataset)
W3 <- prepare_scale_genotypes(plink_data$X, "Asian", new_dataset)
W4 <- prepare_scale_genotypes(plink_data$X, "Mixed", new_dataset)
W5 <- prepare_scale_genotypes(plink_data$X, "Chinese", new_dataset)

# Combine all W matrices vertically, skipping NULL values
W_all <- do.call(rbind, Filter(Negate(is.null), list(W1, W2, W3, W4, W5)))

# Compute the genomic relationship matrix (GRM) using tcrossprod
G <- tcrossprod(W_all) / ncol(W_all)  # ncol now gives the number of SNPs, as we transposed earlier

# Compute eigenvalues and eigenvectors
grm_by_struct_eigen <- eigen(G)

# Assuming row names are already set as individual IDs from 'G'
rownames(grm_by_struct_eigen$vectors) <- rownames(G)

# Save the eigenvalues and eigenvectors to an RDS file
saveRDS(grm_by_struct_eigen, file = opt$output)

