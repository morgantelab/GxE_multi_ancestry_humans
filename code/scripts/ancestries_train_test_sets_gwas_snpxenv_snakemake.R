### Creating training sets ###
rm(list=ls()); gc()
set.seed(123)

# Load necessary libraries
library(optparse)
library(readr)

# Define command-line options
option_list <- list(
  make_option(c("-a", "--ancestry_dir"), type = "character", default = NULL, 
              help = "Directory containing ancestry ID files", metavar = "character"),
  make_option(c("-o", "--output_dir"), type = "character", default = NULL, 
              help = "Directory to save PLINK keep files", metavar = "character"),
  make_option(c("-e", "--exclude_ancestry"), type = "character", default = NULL, 
              help = "Ancestry to exclude for training", metavar = "character")
)

# Parse arguments
opt <- parse_args(OptionParser(option_list = option_list))

# Assign parsed values
ancestry_dir <- opt$ancestry_dir
output_dir <- opt$output_dir
exclude_anc <- opt$exclude_ancestry

# Ensure required arguments are provided
if (is.null(ancestry_dir) | is.null(output_dir) | is.null(exclude_anc)) {
  stop("Both --ancestry_dir, --output_dir, and --exclude_ancestry must be specified.")
}

# Define ancestries
ancestries <- c("asian", "white", "mixed", "black", "chinese")

# Check if the provided ancestry is valid
if (!(exclude_anc %in% ancestries)) {
  stop(paste("Invalid ancestry:", exclude_anc, ". Choose from:", paste(ancestries, collapse=", ")))
}

# Load ancestry ID files into a named list
ancestry_ids <- lapply(ancestries, function(anc) {
  file_path <- file.path(ancestry_dir, paste0(anc, "_ids.txt"))
  read.table(file_path, header = FALSE, stringsAsFactors = FALSE, col.names = c("FID", "IID"))
})
names(ancestry_ids) <- ancestries  # Assign ancestry names to list

# Define the test set: Individuals belonging to the excluded ancestry
test_set <- ancestry_ids[[exclude_anc]]

# Define the training set: All other ancestries combined
train_set <- do.call(rbind, ancestry_ids[names(ancestry_ids) != exclude_anc])

# Ensure output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save files in PLINK format (FID IID)
train_file <- file.path(output_dir, paste0("Train_exclude_", exclude_anc, ".txt"))
test_file <- file.path(output_dir, paste0("Test_", exclude_anc, ".txt"))

write.table(train_set, train_file, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(test_set, test_file, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

cat(paste0("Created Train_exclude_", exclude_anc, ".txt and Test_", exclude_anc, ".txt in ", output_dir, "\n"))
