rm(list=ls()); gc()
set.seed(1123)

library(data.table)

setwd("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/gwas_snp_env")
gwas_df <- fread("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/gwas_snp_env/gxe_gwas_sleep_d_full_dataset.PP0s_interaction_terms.csv")

# Conservative threshold
elbow_logp <- 5.1
p_cutoff <- 10^(-elbow_logp)
subset_conservative <- gwas_df[gwas_df$P < p_cutoff, ]

# Save subset_conservative to CSV
write.csv(subset_conservative, file = "gxe_gwas_sleep_d_PP_subset_conservative_snps_full_datatset_covar.csv", row.names = FALSE)
