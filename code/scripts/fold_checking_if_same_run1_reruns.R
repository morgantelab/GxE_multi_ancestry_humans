### CHecking if two folds are same ###

# Load necessary library
library(readr)

# Define directories
dir1 <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model"
dir2 <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/run_1"

# Number of folds
num_folds <- 5

# Compare folds
for (i in 1:num_folds) {
  file1 <- file.path(dir1, paste0("Fold_", i, ".rds"))
  file2 <- file.path(dir2, paste0("Fold_", i, ".rds"))
  
  # Read files
  fold1 <- readRDS(file1)
  fold2 <- readRDS(file2)
  
  # Check if they are identical
  if (identical(fold1, fold2)) {
    cat(paste0("Fold ", i, ": Identical\n"))
  } else {
    cat(paste0("Fold ", i, ": DIFFERENT\n"))
  }
}
