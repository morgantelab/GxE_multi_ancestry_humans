# Load necessary library
library(dplyr)

load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/data2_20240930.RData")

# Assuming your dataset is named 'df'
dtt <- dtt %>%
  mutate(
    ethn1_consolidated = case_when(
      ethn1_white == 1 ~ "White",
      ethn1_black == 1 ~ "Black",
      ethn1_chinese == 1 ~ "Chinese",
      ethn1_asian == 1 ~ "Asian",
      ethn1_mixed == 1 ~ "Mixed",
      TRUE ~ NA_character_  # if no value matches, set it as NA
    )
  )

save(dtt, file="/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/data2_20240930.RData")
