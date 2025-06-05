# Load libraries
library(optparse)

# Define options
option_list <- list(
  make_option(c("--input_rds"), type = "character", help = "Path to input RDS file"),
  make_option(c("--output_dir"), type = "character", help = "Directory to save output fold RDS files")
)

# Parse arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Load data
data <- readRDS(opt$input_rds)

# Set seed for reproducibility
set.seed(1123)

# Assign folds
n <- nrow(data)
data$Fold <- sample(rep(1:5, length.out = n))
folds <- split(data, data$Fold)

# Save each fold
for (i in 1:5) {
  saveRDS(folds[[i]], file = file.path(opt$output_dir, paste0("Fold_", i, ".rds")))
}
