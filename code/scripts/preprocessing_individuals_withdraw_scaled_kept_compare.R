rm(list = ls()); gc()
set.seed(1123)

library(data.table)

# --- Load withdrawn IDs ---
withdrawn <- fread("/data2/morgante_lab/data/ukbiobank/ind_to_remove/withdraw62347_281_20240209.txt", header = FALSE)$V1

# --- Load missing genotype IDs ---
missing_ids <- fread("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/missing_ids.txt", header = FALSE)$V1

# --- Check overlap: missing_ids ∩ withdrawn ---
overlap_missing <- intersect(withdrawn, missing_ids)
cat(length(overlap_missing), "IDs are in both missing genotype and withdrawn list.\n")
if (length(overlap_missing) > 0) {
  fwrite(data.table(ID = overlap_missing),
         "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/missing_ids_in_withdrawn.txt",
         col.names = FALSE)
  cat("Written to: missing_ids_in_withdrawn.txt\n")
}

# --- Load scaled dataset and extract IDs ---
load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/dec_rerun/scaled_dataset_20250106.Rdata")
scaled_ids <- as.character(dataset$ID)

# --- Check overlap: scaled dataset ∩ withdrawn ---
overlap_scaled <- intersect(withdrawn, scaled_ids)
cat(length(overlap_scaled), "IDs in scaled dataset are also in the withdrawn list.\n")
if (length(overlap_scaled) > 0) {
  fwrite(data.table(ID = overlap_scaled),
         "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/scaled_dataset_withdrawn_overlap.txt",
         col.names = FALSE)
  cat("Written to: scaled_dataset_withdrawn_overlap.txt\n")
}
