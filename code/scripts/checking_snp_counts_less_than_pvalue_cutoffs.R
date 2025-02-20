rm(list=ls()); gc()
set.seed(1123)

# Load necessary libraries
library(data.table)

# Define list of known interaction terms
interaction_terms <- c("Townsend", "act0_d", "TVtime", "sleep_d", "smoking_now", "veg_cook", "fish_oily",
                       "fish_lean", "meat_proc", "poultry", "beef", "lamb", "pork", "cheese", "salt", "tea",
                       "alc1", "waist", "getup", "coffee", "smoked_past", "BFP", "sleep_dev")

# Define the input directory
plink_results_dir <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/gwas_snp_env/"

# Get all PLINK files in the directory
result_files <- list.files(path = plink_results_dir, pattern = "*.glm.linear$", full.names = TRUE)

# Initialize summary table
summary_results <- data.table(File = character(), Total_SNPs = integer(),
                              P_lt_1e5 = integer(), P_lt_1e6 = integer(), P_lt_1e3 = integer(),
                              P_gt_1e5 = integer(), P_gt_1e6 = integer(), P_gt_1e3 = integer(),
                              Check_Total = logical())

# Process each interaction term separately
for (interaction_term in interaction_terms) {
  print(paste("🔹 Processing interaction:", interaction_term))
  
  # Filter files that contain the exact interaction term
  selected_files <- result_files[grepl(paste0("_", interaction_term, "_"), result_files)]
  
  # Skip if no matching files found
  if (length(selected_files) == 0) {
    print(paste("⚠️ No files found for interaction:", interaction_term))
    next
  }
  
  # Process each selected file for this environment
  for (plink_results_file in selected_files) {
    file_name <- basename(plink_results_file)
    print(paste("✅ Processing file:", file_name, "for interaction term:", interaction_term))
    
    # Load PLINK results
    results <- fread(plink_results_file)
    
    # Ensure the required columns exist
    if (!("TEST" %in% colnames(results)) || !("P" %in% colnames(results))) {
      print(paste("⚠️ Skipping", file_name, "- required columns (TEST, P) not found"))
      next
    }
    
    # Ensure P-values are read correctly as numeric
    results[, P := as.numeric(P)]
    
    # Remove NA values in P (if any)
    results <- results[!is.na(P)]
    
    # Construct the correct filter for interaction term
    interaction_test <- paste0("ADDx", interaction_term)
    
    # Subset for interaction term
    results_filtered <- results[TEST == interaction_test]
    
    # Count total SNPs
    total_snps <- nrow(results_filtered)
    
    # Count SNPs by P-value thresholds
    p_lt_1e5 <- nrow(results_filtered[P < 1e-5])
    p_lt_1e6 <- nrow(results_filtered[P < 1e-6])
    p_lt_1e3 <- nrow(results_filtered[P < 1e-3])
    
    p_gt_1e5 <- nrow(results_filtered[P >= 1e-5])
    p_gt_1e6 <- nrow(results_filtered[P >= 1e-6])
    p_gt_1e3 <- nrow(results_filtered[P >= 1e-3])
    
    # Check if counts add up
    check_total <- (p_lt_1e5 + p_gt_1e5 == total_snps) &
      (p_lt_1e6 + p_gt_1e6 == total_snps) &
      (p_lt_1e3 + p_gt_1e3 == total_snps)
    
    # Append results to summary
    summary_results <- rbind(summary_results, data.table(File = file_name, Total_SNPs = total_snps,
                                                         P_lt_1e5 = p_lt_1e5, P_lt_1e6 = p_lt_1e6, P_lt_1e3 = p_lt_1e3,
                                                         P_gt_1e5 = p_gt_1e5, P_gt_1e6 = p_gt_1e6, P_gt_1e3 = p_gt_1e3,
                                                         Check_Total = check_total))
  }
}

# Save summary results
summary_file <- file.path(plink_results_dir, "summary_pvalue_counts.txt")
fwrite(summary_results, summary_file, sep = "\t", quote = FALSE)
print(paste("📊 Summary of SNP counts saved to", summary_file))

print("✅ Processing complete.")
