### Creating training sets ###
rm(list=ls()); gc()
set.seed(123)

# Load necessary library
library(readr)

# Define directories
model_dir <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model"
output_dir <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/gwas_snp_env"

# Number of folds
num_folds <- 5

# Create training and test sets
for (i in 1:num_folds) {
  # Read all folds
  folds <- lapply(1:num_folds, function(j) readRDS(file.path(model_dir, paste0("Fold_", j, ".rds"))))
  
  # Create training and test sets
  test_set <- folds[[i]]  # Leave this fold out for testing
  train_set <- do.call(rbind, folds[-i])  # Combine remaining folds for training
  
  # Extract IDs and create PLINK FID IID format
  test_keep <- data.frame(FID = test_set$ID, IID = test_set$ID)
  train_keep <- data.frame(FID = train_set$ID, IID = train_set$ID)
  
  # Save as PLINK-readable keep files
  write.table(train_keep, file.path(output_dir, paste0("Train_Set_", i, ".txt")), 
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  
  write.table(test_keep, file.path(output_dir, paste0("Test_Set_", i, ".txt")), 
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  
  cat(paste0("Created PLINK keep files for Train_Set_", i, " and Test_Set_", i, " in gwas_snp_env\n"))
}
