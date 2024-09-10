# Load necessary libraries
library(dplyr)
library(optparse)

# Define command-line options
option_list <- list(
  make_option(c("-d", "--dir"), type = "character", default = NULL, 
              help = "path to the working directory", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL, 
              help = "path to the output file", metavar = "character")
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if the required arguments are provided
if (is.null(opt$dir) | is.null(opt$output)) {
  print_help(opt_parser)
  stop("Both the working directory and output file path must be provided", call. = FALSE)
}

# Set the working directory dynamically
setwd(opt$dir)

# Define the ancestries and chromosomes
ancestries <- c("asian", "mixed", "black", "chinese", "white")
chromosomes <- 1:22

# Initialize an empty list to store merged SNPs for each ancestry
merged_snps_by_ancestry <- list()

# Read and merge SNP lists for each ancestry
for (ancestry in ancestries) {
  all_snps <- c()  # Initialize empty vector for each ancestry
  
  for (chr in chromosomes) {
    file_name <- paste0(ancestry, "_filtered_all_chr", chr, ".snplist")
    snp_list <- read.table(file_name, header = FALSE, stringsAsFactors = FALSE)
    snp_ids <- as.character(snp_list$V1)  # Convert to character to avoid any factor conversion issues
    all_snps <- unique(c(all_snps, snp_ids))  # Merge SNPs, ensuring uniqueness
  }
  
  merged_snps_by_ancestry[[ancestry]] <- all_snps
  
  # Print the total number of SNPs for each ancestry
  cat("Total number of SNPs for", ancestry, ":", length(all_snps), "\n")
}

# Initialize the list with SNPs from the first ancestry
common_snps <- merged_snps_by_ancestry[[1]]

# Compare across ancestries
for (i in 2:length(merged_snps_by_ancestry)) {
  common_snps <- intersect(common_snps, merged_snps_by_ancestry[[i]])
  
  # Print intermediate result
  cat("Number of common SNPs after intersecting with", ancestries[i], ":", length(common_snps), "\n")
  
  # Stop if no common SNPs are found
  if (length(common_snps) == 0) {
    cat("No common SNPs found after intersecting with", ancestries[i], "\n")
    break
  }
}

# Count the number of common SNPs
num_common_snps <- length(common_snps)

# Output the final number of common SNPs
cat("Final number of common SNPs across all ancestries:", num_common_snps, "\n")

# Save the common SNPs to a file in the specified directory
writeLines(common_snps, opt$output)
