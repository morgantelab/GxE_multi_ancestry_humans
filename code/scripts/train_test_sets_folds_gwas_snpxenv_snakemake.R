### Creating training sets ###
rm(list=ls()); gc()
set.seed(123)

# Load necessary libraries
library(optparse)
library(readr)

# Define command-line arguments
option_list <- list(
  make_option(c("-m", "--model_dir"), type = "character", default = NULL,
              help = "Directory containing fold .rds files", metavar = "character"),
  make_option(c("-o", "--output_dir"), type = "character", default = NULL,
              help = "Directory to save PLINK txt files", metavar = "character"),
  make_option(c("-n", "--num_folds"), type = "integer", default = 5,
              help = "Number of folds (default: 5)", metavar = "integer")
)

# Parse arguments
opt <- parse_args(OptionParser(option_list = option_list))

# Assign parsed values
model_dir <- opt$model_dir
output_dir <- opt$output_dir
num_folds <- opt$num_folds

# Ensure required arguments are provided
if (is.null(model_dir) | is.null(output_dir)) {
  stop("Both --model_dir and --output_dir must be specified.")
}

# Create training and test sets
for (i in 1:num_folds) {
  # Read all folds
  folds <- lapply(1:num_folds, function(j) readRDS(file.path(model_dir, paste0("Fold_", j, ".rds"))))

  # Create training and test sets
  test_set <- folds[[i]]  # Leave this fold out for testing
  train_set <- do.call(rbind, folds[-i])  # Combine remaining folds for training

  # Extract IDs and create PLINK FID IID format
  test_txt <- data.frame(FID = test_set$ID, IID = test_set$ID)
  train_txt <- data.frame(FID = train_set$ID, IID = train_set$ID)

  # Save as PLINK-readable txt files
  write.table(train_txt, file.path(output_dir, paste0("Train_Set_", i, ".txt")),
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

  write.table(test_txt, file.path(output_dir, paste0("Test_Set_", i, ".txt")),
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

  cat(paste0("Created PLINK txt files for Train_Set_", i, " and Test_Set_", i, " in ", output_dir, "\n"))
}
