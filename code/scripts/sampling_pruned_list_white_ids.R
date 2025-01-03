setwd("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr")

# Load necessary library
set.seed(123) # for reproducibility

# Read the file
id_data <- read.table("pruned_grm_white.rel.id", header = FALSE)

# Randomly sample 20,000 individuals
sampled_ids <- id_data[sample(nrow(id_data), 15000), ]

# Write the sampled IDs to a new file without row or column names
write.table(sampled_ids, "sampled_white_ids.rel.id", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

# Print a message to indicate completion
cat("Sampled 15,000 individuals and saved to sampled_ids.rel.id\n")
