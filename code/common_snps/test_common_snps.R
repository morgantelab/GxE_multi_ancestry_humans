# Load necessary libraries
library(dplyr)

# Define the directory to save sample SNP list files
dir.create("test_snplists", showWarnings = FALSE)

# Sample SNP list data
sample_snps <- list(
  "asian_filtered_all_chr1.snplist" = c("rs1234", "rs2345", "rs3456", "rs4567", "rs5678"),
  "asian_filtered_all_chr2.snplist" = c("rs2345", "rs3456", "rs4567", "rs6789", "rs7890"),
  "mixed_filtered_all_chr1.snplist" = c("rs1234", "rs2345", "rs3456", "rs5678", "rs6789"),
  "mixed_filtered_all_chr2.snplist" = c("rs2345", "rs3456", "rs4567", "rs5678", "rs7890"),
  "black_filtered_all_chr1.snplist" = c("rs1234", "rs3456", "rs4567", "rs5678", "rs7890"),
  "black_filtered_all_chr2.snplist" = c("rs2345", "rs3456", "rs4567", "rs6789", "rs7890"),
  "chinese_filtered_all_chr1.snplist" = c("rs1234", "rs3456", "rs4567", "rs6789", "rs7890"),
  "chinese_filtered_all_chr2.snplist" = c("rs2345", "rs3456", "rs4567", "rs5678", "rs7890")
)

# Write sample SNP list files
for (file_name in names(sample_snps)) {
  writeLines(sample_snps[[file_name]], paste0("test_snplists/", file_name))
}

# Define the ancestries and chromosomes
ancestries <- c("asian", "mixed", "black", "chinese")
chromosomes <- 1:2  # Simplified for testing

# Initialize an empty list to store SNPs
snp_lists <- list()

# Read SNP lists into R
for (ancestry in ancestries) {
  for (chr in chromosomes) {
    file_name <- paste0("test_snplists/", ancestry, "_filtered_all_chr", chr, ".snplist")
    snp_list <- read.table(file_name, header = FALSE, stringsAsFactors = FALSE)
    snp_ids <- as.character(snp_list$V1)  # Convert to character to avoid any factor conversion issues
    snp_lists[[paste0(ancestry, "_chr", chr)]] <- snp_ids
  }
}

# Manual pairwise intersection
common_snps <- snp_lists[[1]]

for (i in 2:length(snp_lists)) {
  common_snps <- intersect(common_snps, snp_lists[[i]])
  
  # Print intermediate result
  cat("Number of common SNPs after intersecting list", i, ":", length(common_snps), "\n")
  
  # Stop if no common SNPs are found
  if (length(common_snps) == 0) {
    cat("No common SNPs found after intersecting list", i, "\n")
    break
  }
}

# Count the number of common SNPs
num_common_snps <- length(common_snps)

# Output the number of common SNPs
cat("Final number of common SNPs across all lists:", num_common_snps, "\n")
