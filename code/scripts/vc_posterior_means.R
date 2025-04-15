#clean enviroment
rm(list = ls())

library(dplyr)
library(tidyr)

setwd("/scratch3/kgoda/ukbiobank_files/tmp/snakemake_runs")

# Define the list of files to process
file_paths <- list(
  # DP = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_DP_pcrelate_just_X1.csv",
  # SP = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_SP_pcrelate_just_X1.csv",
  # PP = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_PP_pcrelate_just_X1.csv"
  # DP = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_DP_pcrelate_just_X1_X2.csv",
  # SP = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_SP_pcrelate_just_X1_X2.csv",
  # # PP = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_PP_pcrelate_just_X1_X2.csv",
  # DP = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_DP_pcrelate_run_G_E_pcrelate_pcs_plink.csv",
  # SP = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_SP_pcrelate_run_G_E_pcrelate_pcs_plink.csv",
  # PP = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_PP_pcrelate_run_G_E_pcrelate_pcs_plink.csv",
  # DP = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/run_1/varabs_DP_pcrelate_pcs_plink_G_run_X_separated.csv",
  # SP = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/run_1/varabs_SP_pcrelate_pcs_plink_G_run_X_separated.csv",
  # PP = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/run_1/varabs_PP_pcrelate_pcs_plink_G_run_X_separated.csv",
  # DP = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/run_1/varabs_DP_pcrelate_run_GE_pcrelate_pcs_plink.csv",
  # SP = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/run_1/varabs_SP_pcrelate_run_GE_pcrelate_pcs_plink.csv",
  # PP = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/run_1/varabs_PP_pcrelate_run_GE_pcrelate_pcs_plink.csv",
  DP = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_DP_pcrelate_X1_X2_E.csv",
  SP = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_SP_pcrelate_X1_X2_E.csv",
  PP = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_PP_pcrelate_X1_X2_E.csv"


)

# Initialize an empty dataframe for storing summary
summary_df <- data.frame()

# Process each file
for (condition in names(file_paths)) {
  # Read the data
  data <- read.csv(file_paths[[condition]], check.names = FALSE)

  # Select only columns V_X1 and V_X2 directly
  data_selected <- data %>% select(any_of(c("V_X1", "V_X2", "V_G", "V_E", "V_GxE")))

  # Calculate means for selected columns using correct syntax
  means_df <- data_selected %>% summarise(across(everything(), \(x) mean(x, na.rm = TRUE))) %>%
    pivot_longer(everything(), names_to = "colname", values_to = "mean")

  # Add filename (condition) to the summary
  means_df <- means_df %>% mutate(filename = condition) %>% select(filename, everything())

  # Append to summary dataframe
  summary_df <- bind_rows(summary_df, means_df)
}

# Define output path
output_file <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/combined_column_means.csv"

# Save summary to CSV
write.csv(summary_df, file = output_file, row.names = FALSE)

cat("Combined column means saved to:", output_file, "\n")
