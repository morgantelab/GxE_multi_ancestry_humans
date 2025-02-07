rm(list=ls()); gc()
set.seed(1123)
library(data.table)

# Define directory where PLINK results are stored
plink_results_dir <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/gwas_snp_env/"

# List of known interaction terms
interaction_terms <- c("Townsend", "act0_d", "TVtime", "sleep_d", "smoking_now", "veg_cook", "fish_oily", 
                       "fish_lean", "meat_proc", "poultry", "beef", "lamb", "pork", "cheese", "salt", "tea", 
                       "alc1", "waist", "getup", "coffee", "smoked_past", "BFP", "sleep_dev")

# Get all PLINK files in the directory
result_files <- list.files(path=plink_results_dir, pattern="*.glm.linear$", full.names=TRUE)

# Dictionary to store results for all interaction terms
ks_results_list <- list()

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
    print(paste("✅ Running KS test for:", file_name, "using interaction term:", interaction_term))
    
    # Load PLINK results
    results <- fread(plink_results_file)
    
    # Define expected column names from PLINK output
    expected_cols <- c("#CHROM", "POS", "ID", "REF", "ALT", "PROVISIONAL_REF?", "A1", "OMITTED", "A1_FREQ", 
                       "TEST", "OBS_CT", "BETA", "SE", "T_STAT", "P", "ERRCODE")
    
    # Define new column names
    new_col_names <- c("CHR", "BP", "SNP", "REF", "ALT", "PROVISIONAL_REF", "A1", "OMITTED", "A1_FREQ", 
                       "TEST", "N", "Beta", "SE", "T", "P", "ERRCODE")
    
    # Identify which columns exist in the PLINK results file
    existing_cols <- intersect(colnames(results), expected_cols)
    
    # Rename only the existing columns
    setnames(results, existing_cols, new_col_names[match(existing_cols, expected_cols)])
    
    # Construct the correct filter for interaction term
    interaction_test <- paste0("ADDx", interaction_term)
    
    # Keep only the SNP × Environment interaction results
    results <- results[TEST == interaction_test]  # Exact filtering
    
    # Ensure P-values are numeric and remove missing values
    results[, P := as.numeric(P)]
    results <- results[!is.na(P)]
    
    # Run KS test comparing p-values to uniform distribution
    if (nrow(results) > 10) {  # Ensure enough SNPs before running KS test
      ks_res <- ks.test(results$P, "punif")
      ks_result <- data.table(File=file_name, Interaction=interaction_term, D=ks_res$statistic, P_value=ks_res$p.value, Num_SNPs=nrow(results))
    } else {
      ks_result <- data.table(File=file_name, Interaction=interaction_term, D=NA, P_value=NA, Num_SNPs=nrow(results))
    }
    
    # Save results for each file
    output_file <- paste0("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/gwas_snp_env/ks_test_result_", interaction_term, "_", file_name, ".txt")
    fwrite(ks_result, output_file, sep="\t")
    
    # Store results in list
    ks_results_list[[paste(interaction_term, file_name, sep="_")]] <- ks_result
  }
}

# Combine all results into a single data table
ks_results <- rbindlist(ks_results_list)

# Save full KS test results
final_output_file <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/gwas_snp_env/ks_test_results_all_envs.txt"
fwrite(ks_results, final_output_file, sep="\t")

print("🎉 KS test results for all environments saved successfully!")
