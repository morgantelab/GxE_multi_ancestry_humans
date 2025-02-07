rm(list=ls()); gc()
library(data.table)
library(ggplot2)

# Define file paths
numeric_results_file <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/gwas_snp_env/gxe_townsend_all_snps_lm_results_numeric.txt"
factor_results_file <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/gwas_snp_env/gxe_townsend_all_snps_lm_results_factor.txt"

# Load results
numeric_results <- fread(numeric_results_file)
factor_results <- fread(factor_results_file)

# Ensure correct column names
setnames(numeric_results, c("Trait", "SNP", "Environment", "Beta", "SE", "P"))
setnames(factor_results, c("Trait", "SNP", "Environment", "Beta", "SE", "P"))

# Merge results by SNP
merged_results <- merge(numeric_results[, .(SNP, P_numeric = P)], 
                        factor_results[, .(SNP, P_factor = P)], 
                        by="SNP", all=TRUE)

# Convert p-values to -log10 scale
merged_results[, logP_numeric := -log10(P_numeric)]
merged_results[, logP_factor := -log10(P_factor)]

# Scatter plot comparing -log10(p-values)
png("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/plots/scatter_sex_numeric_vs_factor.png", width=800, height=600)
ggplot(merged_results, aes(x=logP_numeric, y=logP_factor)) +
  geom_point(alpha=0.5, color="blue") +
  geom_abline(slope=1, intercept=0, color="red", linetype="dashed") +
  labs(title="-log10(P) Comparison: Sex as Numeric vs Factor",
       x="Sex as Numeric (-log10 P)",
       y="Sex as Factor (-log10 P)") +
  theme_minimal()
dev.off()

print("Scatter plot saved successfully!")
