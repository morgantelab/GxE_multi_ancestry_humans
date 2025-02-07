### Creating training sets ###
rm(list=ls()); gc()
set.seed(123)

# Load necessary library
library(readr)

# Define directories
ancestry_dir <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model"
output_dir <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/gwas_snp_env"

# List of ancestries
ancestries <- c("asian", "white", "mixed", "black", "chinese")

# Load ancestry ID files into a named list
ancestry_ids <- lapply(ancestries, function(anc) {
  file_path <- file.path(ancestry_dir, paste0(anc, "_ids.txt"))
  read.table(file_path, header = FALSE, stringsAsFactors = FALSE, col.names = c("FID", "IID"))
})
names(ancestry_ids) <- ancestries  # Assign ancestry names to list

# Create training and test sets excluding one ancestry at a time
for (exclude_anc in ancestries) {
  # Define the test set: Individuals belonging to the excluded ancestry
  test_set <- ancestry_ids[[exclude_anc]]
  
  # Define the training set: All other ancestries combined
  train_set <- do.call(rbind, ancestry_ids[names(ancestry_ids) != exclude_anc])
  
  # Save files in PLINK format (FID IID)
  train_file <- file.path(output_dir, paste0("Train_exclude_", exclude_anc, ".txt"))
  test_file <- file.path(output_dir, paste0("Test_", exclude_anc, ".txt"))
  
  write.table(train_set, train_file, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  write.table(test_set, test_file, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  
  cat(paste0("Created Train_exclude_", exclude_anc, ".txt and Test_", exclude_anc, ".txt\n"))
}
