#clean enviroment
rm(list = ls())

set.seed(1123)

# Load libraries
  library(dplyr)
  library(optparse)

# Define command line options
option_list <- list(
  make_option(c("-i", "--input_rds"), type = "character", help = "Input RDS file path"),
  make_option(c("-l", "--id_list"), type = "character", help = "File containing IDs to retain"),
  make_option(c("-o", "--output_rds"), type = "character", help = "Output RDS file path")
)

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Load dataset
dataset <- readRDS(opt$input_rds)

# Load ID list
id_list <- read.table(opt$id_list, header = FALSE)$V1

# Subset dataset
dataset <- dataset %>% filter(ID %in% id_list)

# Scale traits
dataset$DP0s <- scale(dataset$DP0a) * 10 + 100
dataset$PP0s <- scale(dataset$PP0a) * 10 + 100
dataset$SP0s <- scale(dataset$SP0a) * 10 + 100

# Scale covariates
dataset$AOPs <- scale(dataset$AOP)
dataset$AOPss <- dataset$AOPs^2
dataset$AOPsss <- scale(dataset$AOPss)
dataset$Sex_SIs <- scale(dataset$Sex_SI)

# Save output
saveRDS(dataset, file = opt$output_rds)
