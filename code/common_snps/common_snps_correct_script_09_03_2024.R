rm(list=ls())
setwd("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/output/snp_lists/")


# Load necessary libraries
library(dplyr)

# Define the ancestries and chromosomes
ancestries <- c("asian", "mixed", "black", "chinese")
chromosomes <- 1:22

# Initialize an empty list to store merged SNPs for each ancestry
merged_snps_by_ancestry <- list()

# Read and merge SNP lists for each ancestry
for (ancestry in ancestries) {
  all_snps <- c()
  
  for (chr in chromosomes) {
    file_name <- paste0(ancestry, "_filtered_all_chr", chr, ".snplist")
    snp_list <- read.table(file_name, header = FALSE, stringsAsFactors = FALSE)
    snp_ids <- as.character(snp_list$V1)  # Convert to character to avoid any factor conversion issues
    all_snps <- unique(c(all_snps, snp_ids))  # Merge SNPs, ensuring uniqueness
  }
  
  merged_snps_by_ancestry[[ancestry]] <- all_snps
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
