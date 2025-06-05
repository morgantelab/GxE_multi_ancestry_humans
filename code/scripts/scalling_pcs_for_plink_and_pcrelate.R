#clean enviroment
rm(list = ls())

set.seed(1123)

# Load required libraries
library(optparse)

# Define command line options
option_list <- list(
  make_option(c("--pcrelate_rds"), type = "character", help = "Path to PCA result from PCRelate"),
  make_option(c("--plink_rds"), type = "character", help = "Path to PCA result from PLINK"),
  make_option(c("--out_pcrelate"), type = "character", help = "Output file for scaled PCRelate PCs"),
  make_option(c("--out_plink"), type = "character", help = "Output file for scaled PLINK PCs")
)

# Parse command line options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Read RDS files
eigen_results_pcrelate <- readRDS(opt$pcrelate_rds)
eigen_results_plink <- readRDS(opt$plink_rds)

# Scale PCRelate PCs
ids_pcrelate <- rownames(eigen_results_pcrelate$vectors)
pcs_pcrelate <- eigen_results_pcrelate$vectors[, 1:10]
scaled_pcs_pcrelate <- scale(pcs_pcrelate)
pcs_scaled_pcrelate <- data.frame(ID = ids_pcrelate, scaled_pcs_pcrelate)
saveRDS(pcs_scaled_pcrelate, file = opt$out_pcrelate)

# Scale PLINK PCs
ids_plink <- rownames(eigen_results_plink$vectors)
pcs_plink <- eigen_results_plink$vectors[, 1:10]
scaled_pcs_plink <- scale(pcs_plink)
pcs_scaled_plink <- data.frame(ID = ids_plink, scaled_pcs_plink)
saveRDS(pcs_scaled_plink, file = opt$out_plink)
