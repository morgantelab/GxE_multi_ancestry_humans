# Load necessary libraries
library(SNPRelate)
library(GENESIS)
library(GWASTools)
library(BiocParallel)
library(optparse)

# Define command-line options
option_list <- list(
  make_option(c("-d", "--dir"), type = "character", default = NULL, 
              help = "path to the working directory", metavar = "character"),
  make_option(c("-g", "--gds"), type = "character", default = NULL, 
              help = "path to the GDS file", metavar = "character"),
  make_option(c("-k", "--kinship"), type = "character", default = NULL, 
              help = "path to save the KING kinship matrix", metavar = "character"),
  make_option(c("-p", "--pcair"), type = "character", default = NULL, 
              help = "path to save the PC-AiR results", metavar = "character"),
  make_option(c("-r", "--pcrelate"), type = "character", default = NULL, 
              help = "path to save the PC-Relate results", metavar = "character"),
  make_option(c("-m", "--grm"), type = "character", default = NULL, 
              help = "path to save the GRM matrix", metavar = "character"),
  make_option(c("-c", "--cores"), type = "integer", default = 1, 
              help = "number of cores for parallel computation", metavar = "integer")
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if the required arguments are provided
if (is.null(opt$gds) | is.null(opt$kinship) | is.null(opt$pcair) | is.null(opt$pcrelate) | is.null(opt$grm)) {
  print_help(opt_parser)
  stop("All file paths must be provided.", call. = FALSE)
}

# Set the working directory dynamically
setwd(opt$dir)

# Open the GDS file
gds <- snpgdsOpen(opt$gds)

# Calculate KING-robust kinship coefficients
king_kinship <- snpgdsIBDKING(
  gds,
  sample.id = NULL,           # Include all samples
  snp.id = NULL,              # Include all SNPs
  autosome.only = TRUE,       # Use only autosomal SNPs
  remove.monosnp = TRUE,      # Do not remove monomorphic SNPs
  maf = NaN,                  # No filtering by MAF
  missing.rate = NaN,         # No filtering by missing rate
  type = "KING-robust",       # Use KING-robust kinship estimation
  family.id = NULL,           # No family IDs
  num.thread = opt$cores,     # Number of threads for parallel computation
  verbose = TRUE              # Print detailed information
)

# Change up the KING matrix for further analysis
king_kinship_mat <- king_kinship$kinship
colnames(king_kinship_mat) <- rownames(king_kinship_mat) <- king_kinship$sample.id

# Save the KING kinship matrix
save(king_kinship_mat, file = opt$kinship)

# Close the GDS file after analysis
snpgdsClose(gds)

# Load genotype data for PCA analysis
diversity_geno <- GdsGenotypeReader(filename = opt$gds)
diversity_genoData <- GenotypeData(diversity_geno)

# Run PC-AiR
pcair_results <- pcair(
  diversity_genoData,
  kinobj = king_kinship_mat,  # KING kinship matrix
  divobj = king_kinship_mat,  # Use the same matrix for both kin and div
  num.cores = opt$cores,      # Use the number of cores provided
  verbose = TRUE
)

# Save the PC-AiR results
saveRDS(pcair_results, file = opt$pcair)

# Load PC-AiR results for PC-Relate
diversity_genoData <- GenotypeBlockIterator(diversity_genoData)

# Run PC-Relate with 5 principal components
mypcrelate_5pcs <- pcrelate(
  diversity_genoData, 
  pcs = pcair_results$vectors[, 1:5],  # Use first 5 PCs
  ibd.probs = FALSE, 
  training.set = pcair_results$unrels, 
  BPPARAM = BiocParallel::MulticoreParam(workers = opt$cores)
)

# Save the PC-Relate results
saveRDS(mypcrelate_5pcs, file = opt$pcrelate)

# Extract individual IDs
iids <- as.character(getScanID(diversity_genoData))

# Create GRM matrix from PC-Relate results
grm_matrix_pcrelate_5pcs <- pcrelateToMatrix(mypcrelate_5pcs, sample.include = iids, thresh = NULL)

# Save the GRM matrix
save(grm_matrix_pcrelate_5pcs, file = opt$grm)

# Close genotype data objects
close(diversity_geno)
print("Analysis completed.")

