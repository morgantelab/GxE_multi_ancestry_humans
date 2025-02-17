rm(list=ls()); gc()
set.seed(1123)

# Load required libraries
library(data.table)
library(optparse)

# Define command-line options
option_list <- list(
  make_option(c("-d", "--dir"), type = "character", default = NULL,
              help = "Path to the working directory", metavar = "character"),
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Path to the input GWAS file (.glm.linear)", metavar = "character"),
  make_option(c("-e", "--env"), type = "character", default = NULL,
              help = "environment covariate", metavar = "character"),
  make_option(c("-f", "--fold"), type = "character", default = NULL,
              help = "fold", metavar = "character"),
  make_option(c("-t", "--trait"), type = "character", default = NULL,
              help = "trait", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output directory for saving filtered results", metavar = "character")
)

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check required arguments
if (is.null(opt$dir) | is.null(opt$input) | is.null(opt$output) | is.null(opt$trait) | is.null(opt$env) | is.null(opt$fold)) {
  print_help(opt_parser)
  stop("All file paths must be provided.", call. = FALSE)
}

# Set working directory
setwd(opt$dir)

# Read input GWAS file
cat("Processing file:", opt$input, "\n")
df <- fread(opt$input, header = TRUE)

# Identify the P-value column dynamically
p_col <- grep("P$", colnames(df), value = TRUE)
if (length(p_col) == 0) {
  stop("No P-value column found in", opt$input)
}

# Extract environmental variable from filename
file_basename <- basename(opt$input)
env_var <- sub("gxe_([^_]+)_.*", "\\1", file_basename)  # Extract env name after 'gxe_'

# Filter for interaction terms (ADDx{env}) where P-value < 1e-5
df_filtered <- df[grepl(paste0("ADDx", env_var, "$"), df$TEST) & df[[p_col]] < 1e-5, ]

if (nrow(df_filtered) == 0) {
  cat("No significant interactions found for", env_var, "in", opt$input, "\n")
  quit(save = "no")
}

# Define output filename
output_file <- file.path(opt$output, paste0(file_basename, ".filtered.linear"))

# Save filtered results
fwrite(df_filtered, output_file, sep = "\t")

cat("Filtered results saved to:", output_file, "\n")
