#### Sex and Age added to Emat, Emat creation, Processing Emat, Eigen of Emat. ####

rm(list=ls()); gc()

# Set seed for reproducibility
set.seed(1123)

# Load necessary libraries
library(optparse)

# Define command-line options
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Path to input RData file containing dataset", metavar="character"),
  make_option(c("-d", "--outdir"), type="character", default="output/",
              help="Directory where output files will be saved", metavar="character")
)

# Parse options
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate input file
if (is.null(opt$input)) {
  stop("Error: Input file is required. Use -i or --input to specify the dataset RData file.")
}

# Define output file paths
emat_path <- file.path(opt$outdir, "Emat_with_sex_age.RDS")
eigen_path <- file.path(opt$outdir, "Emat_eigen_with_sex_age.RDS")

# Load dataset
load(opt$input)

# Check if dataset exists
if (!exists("dataset")) {
  stop("Error: 'dataset' object not found in the input file.")
}

# Compute sleep deviation
dataset$sleep_dev <- (dataset$sleep - mean(dataset$sleep))^2

# Define environmental variables, EXCLUDING AOPs and AOPss - ran but not used. output in Emat_with_sex_age_not_used folder in data.
# Evars <- c("Sex_SI", "Townsend", "act0_d", "TVtime", "sleep_d",
#            "smoking_now", "veg_cook", "fish_oily", "fish_lean", "meat_proc",
#            "poultry", "beef", "lamb", "pork", "cheese", "salt", "tea", "alc1",
#            "waist", "getup", "coffee", "smoked_past", "BFP", "sleep_dev")

# Define environmental variables, EXCLUDING AOPs and AOPss
Evars <- c("AOPss", "Sex_SI", "Townsend", "act0_d", "TVtime", "sleep_d",
           "smoking_now", "veg_cook", "fish_oily", "fish_lean", "meat_proc",
           "poultry", "beef", "lamb", "pork", "cheese", "salt", "tea", "alc1",
           "waist", "getup", "coffee", "smoked_past", "BFP", "sleep_dev")

# Check if all Evars exist in dataset
missing_vars <- setdiff(Evars, names(dataset))
if (length(missing_vars) > 0) {
  stop(paste("Error: The following variables are missing in the dataset:", paste(missing_vars, collapse=", ")))
}

# Convert selected columns to matrix and scale
Emat_scaled <- scale(as.matrix(dataset[Evars]))

# # Extract already scaled AOPs and AOPss - ran but not used. output in Emat_with_sex_age_not_used folder in data.
# if (!all(c("AOPs", "AOPss") %in% names(dataset))) {
#   stop("Error: AOPs or AOPss not found in dataset.")
# }
# AOP_mat <- as.matrix(dataset[, c("AOPs", "AOPss")])  # Use as.matrix to ensure proper structure

# Extract already scaled AOPs
if (!"AOPs" %in% names(dataset)) {
  stop("Error: AOPs not found in dataset.")
}
AOP_mat <- as.matrix(dataset[,"AOPs", drop=FALSE])  # Ensure it's a matrix

# Combine scaled Emat with unscaled AOPs and AOPss
Emat <- cbind(AOP_mat, Emat_scaled)

# Assign row names
rownames(Emat) <- dataset$ID

# Compute E matrix
E <- tcrossprod(Emat)
E <- E / mean(diag(E))

# Compute eigenvalues and eigenvectors
E_eigen <- eigen(E)

# Set row names
rownames(E_eigen$vectors) <- rownames(Emat)

# Save outputs
saveRDS(E_eigen, file=eigen_path)
save(Emat, file=emat_path)

cat("Eigenvalues and eigenvectors saved to:", eigen_path, "\n")
cat("Scaled environmental matrix (Emat) saved to:", emat_path, "\n")
