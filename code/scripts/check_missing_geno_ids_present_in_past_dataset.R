##quick check of missing_geno_ids if present in used dataset ##

rm(list = ls()); gc()

library(data.table)


# Load current dataset
load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/dec_rerun/scaled_dataset_20250106.Rdata")
current_ids <- as.character(dataset$ID)

# Load missing IDs
missing_ids <- fread("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/missing_ids.txt", header = FALSE)$V1

# Check overlap
intersect_ids <- intersect(current_ids, missing_ids)
cat(length(intersect_ids), "IDs with missing genotypes are present in your current dataset.\n")

# Optional: write to file if needed
if (length(intersect_ids) > 0) {
  write.table(intersect_ids,
              "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/missing_ids_in_scaled_dataset.txt",
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  cat("Written to: missing_ids_in_scaled_dataset.txt\n")
}