rm(list=ls()); gc()
library(data.table)
library(qqman)
library(ggplot2)
library(genio)  # Use genio to read PLINK files

# Define file paths
results_file <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/gwas_snp_env/gxe_townsend_all_snps_lm_results.txt"
plink_prefix <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/merged_geno_common_snps_selected_indiv"

# Load GWAS results
results <- fread(results_file)

# Load SNP metadata using `genio`
plink_data <- read_plink(plink_prefix)  # Reads .bim, .bed, .fam files
snp_metadata <- plink_data$bim  # Extract SNP metadata

# Ensure SNP metadata has correct column names
setnames(snp_metadata, c("CHR", "SNP", "CM", "BP", "A1", "A2"))

# Merge GWAS results with SNP metadata
results <- merge(results, snp_metadata[, c("CHR", "SNP", "BP")], by="SNP", all.x=TRUE)

# Remove missing values and ensure correct data types
results <- results[!is.na(CHR) & !is.na(BP)]
results[, P := as.numeric(P)]

# Sort results for Manhattan plotting
setorder(results, CHR, BP)

# Convert CHR to numeric and remove non-autosomal values
results[, CHR := as.integer(as.character(CHR))]

# Check if conversion worked
if (any(is.na(results$CHR))) {
  print("Warning: Some CHR values could not be converted.")
  print(table(results$CHR, useNA="always"))  # Print CHR distribution
}

# Remove non-autosomal chromosomes (e.g., X, Y, MT if they exist)
results <- results[CHR %in% 1:22]  # Keep only chromosomes 1-22

# **Manhattan Plot for SNP × Townsend Interaction**
png("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/plots/manhattan_gxe_townsend.png", width=1200, height=600)
manhattan(results, chr="CHR", bp="BP", p="P", snp="SNP",
          main="Manhattan Plot: SNP × Townsend Interaction",
          ylim=c(0, -log10(min(results$P, na.rm=TRUE)) + 1),
          col=c("blue4", "orange3"), suggestiveline=-log10(1e-5), genomewideline=-log10(5e-8))
dev.off()

# **QQ Plot for SNP × Environment Interaction**
png("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/plots/qqplot_gxe_townsend.png", width=800, height=600)
qq(results$P, main="QQ Plot: SNP × Townsend Interaction")
dev.off()

print("Plots saved successfully!")
