
rm(list=ls()); gc()
set.seed(1123)

# Load required libraries
library(genio)

# Load PLINK data using genio
plink_data <- read_plink("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/merged_geno_common_snps_selected_indiv")

# Extract SNP genotype matrix (Missing genotypes are stored as NA)
snp_matrix <- plink_data$X

# Count missing genotypes per individual (rows)
missing_per_snp <- rowSums(is.na(snp_matrix))

# Count missing genotypes per SNP (columns)
missing_per_individual <- colSums(is.na(snp_matrix))

# Check number of individuals (rows) and SNPs (columns)
num_individuals <- nrow(plink_data$fam)
num_snps <- nrow(plink_data$bim)

# Compute the expected number of missing genotypes
total_genotypes <- as.numeric(num_individuals) * as.numeric(num_snps)
observed_non_missing <- sum(!is.na(snp_matrix))  # Count of non-missing genotypes
expected_missing_count <- total_genotypes - observed_non_missing  # Expected missing count

# Print summary
cat("Number of SNPs:", num_snps, "\n")
cat("Number of individuals:", num_individuals, "\n")
cat("Total genotype calls (expected):", total_genotypes, "\n")
cat("Observed non-missing genotype calls:", observed_non_missing, "\n")
cat("Total missing genotypes (expected from PLINK 01 encoding):", expected_missing_count, "\n\n")

# Print summary statistics for missing genotypes
cat("Summary of missing genotypes per individual:\n")
print(summary(missing_per_individual))

cat("\nSummary of missing genotypes per SNP:\n")
print(summary(missing_per_snp))

# Optionally, export missing genotype information
write.csv(missing_per_individual, "missing_per_individual.csv", row.names = TRUE)
write.csv(missing_per_snp, "missing_per_snp.csv", row.names = TRUE)
