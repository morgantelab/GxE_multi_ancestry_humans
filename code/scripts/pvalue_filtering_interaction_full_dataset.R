rm(list=ls()); gc()
set.seed(1123)
library(data.table)

# Define directories
results_dir <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/gwas_snp_env/"
output_dir <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/gwas_snp_env/"

# List of environments for which interaction terms are being analyzed
envs <- c("sex", "age", "age2")
# envs <- c("Townsend", "act0_d", "TVtime", "sleep_d", "smoking_now",
#           "veg_cook", "fish_oily", "fish_lean", "meat_proc", "poultry",
#           "beef", "lamb", "pork", "cheese", "salt", "tea", "alc1",
#           "waist", "getup", "coffee", "smoked_past", "BFP", "sleep_dev")


# Loop through each environment
for (env in envs) {
  print(paste("Processing Environment:", env))

  # Find all files corresponding to this environment
  env_files <- list.files(results_dir, pattern=paste0("S_A_gxe_", env, "_full_dataset.*\\.glm\\.linear$"), full.names=TRUE)

  for (file in env_files) {
    print(paste("Reading file:", file))

    # Read GWAS results
    plink_results <- fread(file)

    # Rename columns if necessary
    setnames(plink_results, c("#CHROM", "POS", "ID", "REF", "ALT", "A1", "TEST", "OBS_CT", "BETA", "SE", "T_STAT", "P", "ERRCODE"),
             c("CHR", "BP", "SNP", "REF", "ALT", "A1", "TEST", "N", "Beta", "SE", "T", "P", "ERRCODE"), skip_absent=TRUE)

    # Extract only the interaction term for this environment
    interaction_term <- paste0("ADDx", env)
    interaction_results <- plink_results[TEST == interaction_term]

    # Convert CHR to numeric and filter autosomal chromosomes
    interaction_results[, CHR := as.integer(as.character(CHR))]
    interaction_results <- interaction_results[CHR %in% 1:22]

    # Ensure P-values are numeric
    interaction_results[, P := as.numeric(P)]

    # Sort the results by p-value (smallest first)
    interaction_results <- interaction_results[order(P)]

    # Generate an output filename based on the input file
    output_filename <- gsub("\\.glm\\.linear$", "_interaction_terms.csv", basename(file))
    output_filepath <- file.path(output_dir, output_filename)

    # Save the sorted interaction terms to a separate file
    if (nrow(interaction_results) > 0) {
      fwrite(interaction_results, output_filepath)
      print(paste("Saved interaction terms to:", output_filepath))
    } else {
      print(paste("No interaction terms found for:", file))
    }
  }
}

print("All interaction term files have been processed and saved.")


