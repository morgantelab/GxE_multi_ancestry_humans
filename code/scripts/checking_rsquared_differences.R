# Load necessary library
library(dplyr)

# Set base directory
base_dir <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/"

# Define file paths
#file1 <- "R2_Corr_Results_ethnicities (2).csv"
#file2 <- "all_model_summary_pval_ethn_corrected.csv"

# Read in the data
#df1 <- read.csv(file1, stringsAsFactors = FALSE)
#df2 <- read.csv(file2, stringsAsFactors = FALSE)

df1 <- R2_Corr_Results_ethnicities_2_
df2 <- all_model_summary_pval_ethn_corrected

# Ensure column names match
colnames(df1) <- trimws(colnames(df1))
colnames(df2) <- trimws(colnames(df2))

# Convert Ancestry column to lowercase for consistent matching
df1$Ancestry <- tolower(df1$Ancestry)
df2$Ancestry <- tolower(df2$Ancestry)

# Filter out rows where Pval is not NA in df2
df2_filtered <- df2 %>% filter(is.na(Pval))

# Merge the dataframes on Model, Ancestry, and Trait
merged_df <- df1 %>%
  inner_join(df2_filtered, by = c("Model", "Ancestry", "Trait"), suffix = c("_file1", "_file2"))

# Compute differences in R_squared
merged_df <- merged_df %>%
  mutate(R_squared_diff = R_squared_file1 - R_squared_file2)

# Select relevant columns
result <- merged_df %>%
  select(Model, Ancestry, Trait, R_squared_file1, R_squared_file2, R_squared_diff)

# Print the results
print(result)

# Optionally, save to CSV
write.csv(result, "R_squared_differences_past_present_plots_also_seen.csv", row.names = FALSE)
