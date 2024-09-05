rm(list=ls())
setwd("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/output/snp_lists/")

# Load necessary libraries
library(dplyr)

# Define the ancestries and chromosomes
ancestries <- c("asian", "mixed", "black", "chinese")
chromosomes <- 1:22

# Initialize an empty vector to store all SNPs
all_snps <- c()

# Read SNP lists into R and combine them
for (ancestry in ancestries) {
  for (chr in chromosomes) {
    file_name <- paste0(ancestry, "_filtered_all_chr", chr, ".snplist")
    snp_list <- read.table(file_name, header = FALSE, stringsAsFactors = FALSE)
    snp_ids <- as.character(snp_list$V1)  # Convert to character to avoid any factor conversion issues
    all_snps <- c(all_snps, snp_ids)  # Combine all SNPs into one vector
  }
}

# Count occurrences of each SNP
snp_table <- table(all_snps)

# Filter SNPs that appear in all lists
num_lists <- length(ancestries) * length(chromosomes)  # Total number of SNP lists
common_snps <- names(snp_table[snp_table == num_lists])

# Count the number of common SNPs
num_common_snps <- length(common_snps)

# Output the number of common SNPs
cat("Number of common SNPs across all lists:", num_common_snps, "\n")


