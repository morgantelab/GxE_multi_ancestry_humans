rm(list=ls()); gc()
set.seed(1123)

# Load required library
library(Matrix)

# Load both Hadamard GRM files
load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/gwas_snp_env/Hadamard_GRM_Fold1_SP_samplerun.RData")
grm1 <- K_fold_grm  # Assuming the object is named Hadamard_GRM

load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/gwas_snp_env/Hadamard_GRM_Fold1_SP.RData")
grm2 <- K_fold_grm  # Assuming the object is named Hadamard_GRM

# Check if both GRMs are identical
if (identical(grm1, grm2)) {
  cat("The two Hadamard GRMs are exactly the same.\n")
} else {
  cat("The two Hadamard GRMs are NOT the same.\n")
  
  # Optional: Check for element-wise differences
  if (all(dim(grm1) == dim(grm2))) {
    diff_matrix <- grm1 - grm2
    max_diff <- max(abs(diff_matrix))
    
    if (max_diff == 0) {
      cat("The matrices have the same dimensions and values, but might differ in attributes.\n")
    } else {
      cat("Maximum absolute difference between elements:", max_diff, "\n")
    }
  } else {
    cat("The matrices have different dimensions.\n")
  }
}
