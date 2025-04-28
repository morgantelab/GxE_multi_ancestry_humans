rm(list=ls()); gc()
set.seed(1123)

library(data.table)

setwd("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/gwas_snp_env")
gwas_df <- fread("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/gwas_snp_env/gxe_gwas_BFP_full_dataset.DP0s_interaction_terms.csv")

# Conservative threshold
subset_conservative <- gwas_df[gwas_df$P < 1e-5, ]

# Save subset_conservative to CSV
write.csv(subset_conservative, file = "gxe_gwas_BFP_DP_subset_conservative_snps_full_datatset_covar.csv", row.names = FALSE)
