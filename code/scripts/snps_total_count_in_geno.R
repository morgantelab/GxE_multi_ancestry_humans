# Clear workspace
rm(list=ls())
gc()

library(genio)

# Specify the path to your bed file (without the extension)
bed_prefix <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/merged_geno_common_snps_selected_indiv"

# Read the .bim file to efficiently obtain SNP information
bim_data <- read_bim(bed_prefix)

# Count the number of SNPs
num_snps <- nrow(bim_data)

# Print the result
cat("Total number of SNPs:", num_snps, "\n")
