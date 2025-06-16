rm(list=ls()); gc()
library(data.table)
library(ggplot2)

# Define file paths
lm_results_file <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/gwas_snp_env/gxe_townsend_all_snps_lm_results_factor.txt"
plink_results_file <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/gwas_snp_env/plink_gwas_snp_1env_test_param_1_16.DP0s.glm.linear"

# Load `lm()` results
lm_results <- fread(lm_results_file)
setnames(lm_results, c("Trait", "SNP", "Environment", "Beta", "SE", "P"))

# Load PLINK results
plink_results <- fread(plink_results_file)

# Rename columns based on actual PLINK output
setnames(plink_results, c("#CHROM", "POS", "ID", "REF", "ALT", "PROVISIONAL_REF?", "A1", "OMITTED", "A1_FREQ", 
                          "TEST", "OBS_CT", "BETA", "SE", "T_STAT", "P", "ERRCODE"),
         c("CHR", "BP", "SNP", "REF", "ALT", "PROVISIONAL_REF", "A1", "OMITTED", "A1_FREQ", 
           "TEST", "N", "Beta", "SE", "T", "P", "ERRCODE"))

# Keep only SNP × Environment interaction term (adjust TEST as needed)
plink_results <- plink_results[TEST == "ADDxTownsend"]

# Merge results by SNP
merged_results <- merge(lm_results[, .(SNP, P_lm = P)], 
                        plink_results[, .(SNP, P_plink = P)], 
                        by="SNP", all=TRUE)

# Convert p-values to -log10 scale
merged_results[, logP_lm := -log10(P_lm)]
merged_results[, logP_plink := -log10(P_plink)]

# Scatter plot comparing -log10(p-values)
png("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/plots/scatter_plink_vs_lm_param.png", width=800, height=600)
ggplot(merged_results, aes(x=logP_lm, y=logP_plink)) +
  geom_point(alpha=0.5, color="blue") +
  geom_abline(slope=1, intercept=0, color="red", linetype="dashed") +
  labs(title="-log10(P) Comparison: PLINK vs. lm()",
       x="lm() -log10(P)",
       y="PLINK -log10(P)") +
  theme_minimal()
dev.off()

print("Scatter plot saved successfully!")
