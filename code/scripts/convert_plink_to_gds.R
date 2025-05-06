set.seed(1123)
# Load the SNPRelate package
library(SNPRelate)
library(optparse)

# Define command-line options
option_list <- list(
  make_option(c("-d", "--dir"), type = "character", default = NULL,
              help = "path to the working directory", metavar = "character"),
  make_option(c("-b", "--bed"), type = "character", default = NULL,
              help = "path to the PLINK .bed file", metavar = "character"),
  make_option(c("-i", "--bim"), type = "character", default = NULL,
              help = "path to the PLINK .bim file", metavar = "character"),
  make_option(c("-f", "--fam"), type = "character", default = NULL,
              help = "path to the PLINK .fam file", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "path to the output GDS file", metavar = "character")
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if the required arguments are provided
if (is.null(opt$dir) | is.null(opt$output) | is.null(opt$bed) | is.null(opt$bim) | is.null(opt$fam)) {
  print_help(opt_parser)
  stop("file path must be provided", call. = FALSE)
}

# Set the working directory dynamically
setwd(opt$dir)

# Specify the paths to your PLINK files
#bed_file <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/apr_run_with_chinese_20240402/2_sampling_then_grm_20240402/7_grm_genesis/selected_combined_geno.bed"
#bim_file <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/apr_run_with_chinese_20240402/2_sampling_then_grm_20240402/7_grm_genesis/selected_combined_geno.bim"
#fam_file <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/apr_run_with_chinese_20240402/2_sampling_then_grm_20240402/7_grm_genesis/selected_combined_geno.fam"

# Specify the output GDS file name
#output_gds_file <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/apr_run_with_chinese_20240402/2_sampling_then_grm_20240402/7_grm_genesis/selected_combined_geno.gds"

# Convert PLINK to GDS
snpgdsBED2GDS(bed.fn = opt$bed, bim.fn = opt$bim, fam.fn = opt$fam, out.gdsfn = opt$output)

# Inform the user the conversion is done
cat("GDS file created at:", opt$output, "\n")
