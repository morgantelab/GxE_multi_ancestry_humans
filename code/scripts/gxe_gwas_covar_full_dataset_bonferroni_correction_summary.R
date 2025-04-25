rm(list=ls()); gc()
set.seed(1123)

library(data.table)
library(stringr)

# Define directories
results_dir <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/gwas_snp_env/"
output_dir <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/gwas_snp_env/"

# List of environments for which interaction terms are being analyzed
#envs <- c("sex", "age", "age2")
envs <- c("Townsend", "act0_d", "TVtime", "sleep_d", "smoking_now",
          "veg_cook", "fish_oily", "fish_lean", "meat_proc", "poultry",
          "beef", "lamb", "pork", "cheese", "salt", "tea", "alc1",
          "waist", "getup", "coffee", "smoked_past", "BFP", "sleep_dev")

# Define parameters
n_snps <- 302688
n_envs <- length(envs)
n_traits <- 3
total_tests <- n_snps
bonf_threshold <- 0.05 / total_tests

cat("Bonferroni threshold:\n", bonf_threshold, "\n\n")

# Store summary
summary_list <- list()

# Loop through each environment
for (env in envs) {
  cat("Processing Environment:", env, "\n")
  
  env_files <- list.files(results_dir, pattern = paste0("gxe_gwas_", env, "_full_dataset.*\\.glm\\.linear$"), full.names = TRUE)
  
  for (file in env_files) {
    trait <- ifelse(str_detect(file, "DP0s"), "DBP",
                    ifelse(str_detect(file, "SP0s"), "SBP",
                           ifelse(str_detect(file, "PP0s"), "PP", NA)))
    
    if (is.na(trait)) next  # Skip if trait not identified
    
    plink_results <- fread(file)
    setnames(plink_results,
             c("#CHROM", "POS", "ID", "REF", "ALT", "A1", "TEST", "OBS_CT", "BETA", "SE", "T_STAT", "P", "ERRCODE"),
             c("CHR", "BP", "SNP", "REF", "ALT", "A1", "TEST", "N", "Beta", "SE", "T", "P", "ERRCODE"),
             skip_absent = TRUE)
    
    interaction_term <- paste0("ADDx", env)
    interaction_results <- plink_results[TEST == interaction_term]
    interaction_results[, P := as.numeric(P)]
    
    n_bonf <- sum(interaction_results$P < bonf_threshold, na.rm = TRUE)
    summary_list[[length(summary_list) + 1]] <- data.table(Environment = env,
                                                           Trait = trait,
                                                           N_Bonferroni = n_bonf)
  }
}

# Combine summary
summary_df <- rbindlist(summary_list)
summary_wide <- dcast(summary_df, Environment ~ Trait, value.var = "N_Bonferroni", fill = 0)

# View result
print(summary_wide)

# Optionally save
fwrite(summary_wide, file.path(output_dir, "gxe_gwas_bonferroni_pass_summary_by_trait_covar_full_dataset.csv"))
