# Clean environment
rm(list = ls())
set.seed(1123)
library(dplyr)
library(tidyr)

setwd("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model")

# Define the list of files to process
file_paths <- list(
  DP_X1 = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_DP_pcrelate_just_X1.csv",
  SP_X1 = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_SP_pcrelate_just_X1.csv",
  PP_X1 = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_PP_pcrelate_just_X1.csv",
  DP_X1_X2 = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_DP_pcrelate_just_X1_X2.csv",
  SP_X1_X2 = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_SP_pcrelate_just_X1_X2.csv",
  PP_X1_X2 = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_PP_pcrelate_just_X1_X2.csv",
  DP_G_E = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_DP_pcrelate_X1_X2_G_E.csv",
  SP_G_E = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_SP_pcrelate_X1_X2_G_E.csv",
  PP_G_E = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_PP_pcrelate_X1_X2_G_E.csv",
  DP_G = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_DP_pcrelate_X1_X2_G.csv",
  SP_G = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_SP_pcrelate_X1_X2_G.csv",
  PP_G = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_PP_pcrelate_X1_X2_G.csv",
  DP_GxE = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_DP_pcrelate_X1_X2_G_E_GE.csv",
  SP_GxE = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_SP_pcrelate_X1_X2_G_E_GE.csv",
  PP_GxE = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_PP_pcrelate_X1_X2_G_E_GE.csv",
  DP_E = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_DP_pcrelate_X1_X2_E.csv",
  SP_E = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_SP_pcrelate_X1_X2_E.csv",
  PP_E = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_PP_pcrelate_X1_X2_E.csv"
)

# Initialize an empty dataframe for storing summary
summary_df <- data.frame()

# Process each file
for (condition in names(file_paths)) {
  # Read the data
  data <- read.csv(file_paths[[condition]], check.names = FALSE)

  # Select relevant columns
  data_selected <- data %>% select(any_of(c("V_X1", "V_X2", "V_E", "V_G", "V_GE")))

  # Reshape to long format for easy summary
  data_long <- data_selected %>%
    pivot_longer(cols = everything(), names_to = "colname", values_to = "value")

  # Calculate mean and 95% CI for each component
  stats <- data_long %>%
    group_by(colname) %>%
    summarise(
      mean = mean(value, na.rm = TRUE),
      ci_lower = mean(value, na.rm = TRUE) - qt(0.975, sum(!is.na(value)) - 1) * sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))),
      ci_upper = mean(value, na.rm = TRUE) + qt(0.975, sum(!is.na(value)) - 1) * sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))),
      .groups = "drop"
    ) %>%
    mutate(filename = condition) %>%
    select(filename, everything())

  # Append to final summary
  summary_df <- bind_rows(summary_df, stats)
}

# Save summary to file
output_file <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/vc_means_all_models_with_ci.csv"
write.csv(summary_df, file = output_file, row.names = FALSE)

cat("Combined column means with confidence intervals saved to:", output_file, "\n")
