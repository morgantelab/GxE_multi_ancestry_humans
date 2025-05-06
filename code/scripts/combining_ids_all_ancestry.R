set.seed(1123)

setwd("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr")

# Read each ID list
sampled_white_ids <- read.table("sampled_white_ids.rel.id", header = FALSE)
pruned_grm_asian <- read.table("pruned_grm_asian.rel.id", header = FALSE)
pruned_grm_black <- read.table("pruned_grm_black.rel.id", header = FALSE)
pruned_grm_chinese <- read.table("pruned_grm_chinese.rel.id", header = FALSE)
pruned_grm_mixed <- read.table("pruned_grm_mixed.rel.id", header = FALSE)

# Merge all the ID lists into one data frame
merged_ids <- rbind(sampled_white_ids, pruned_grm_asian, pruned_grm_black, pruned_grm_chinese, pruned_grm_mixed)

# Remove duplicates if any (optional, if needed)
merged_ids <- unique(merged_ids)

# Write the merged list to a new file without row or column names
write.table(merged_ids, "merged_ids.rel.id", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

# Print a message to indicate completion
cat("Merged ID lists and saved to merged_ids.rel.id\n")
