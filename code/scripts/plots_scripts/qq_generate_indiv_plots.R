rm(list=ls()); gc()

set.seed(1123)

library(data.table)
library(ggplot2)
library(qqman)

# Define directories
results_dir <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/gwas_snp_env/"
output_dir <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/plots/"

#List of 23 environments
envs <- c("Townsend", "act0_d", "TVtime", "sleep_d", "smoking_now",
          "veg_cook", "fish_oily", "fish_lean", "meat_proc", "poultry",
          "beef", "lamb", "pork", "cheese", "salt", "tea", "alc1",
          "waist", "getup", "coffee", "smoked_past", "BFP", "sleep_dev")

#envs <- c("sex", "age", "age2")

# Function to generate a QQ plot
generate_qq_plot <- function(p_values, title, output_file) {
  p_values <- as.numeric(p_values)  # Ensure numeric conversion
  p_values <- p_values[!is.na(p_values) & is.finite(p_values) & p_values > 0]  # Filter valid values

  if (length(p_values) == 0) {
    print(paste("Skipping", title, "- No valid p-values"))
    return(NULL)  # Skip plot if no valid data
  }

  png(output_file, width=800, height=600)
  qq(p_values, main=title)
  dev.off()
}

# Loop through each environment
for (env in envs) {
  print(paste("Processing Environment:", env))

  # Find all files corresponding to this environment
  env_files <- list.files(results_dir, pattern=paste0("gxe_gwas_", env, "_full_dataset.*\\.glm\\.linear$"), full.names=TRUE)

  for (file in env_files) {
    print(paste("Reading file:", file))

    # Read GWAS results
    plink_results <- fread(file)

    # Rename columns based on actual PLINK output
    setnames(plink_results, c("#CHROM", "POS", "ID", "REF", "ALT", "PROVISIONAL_REF?", "A1", "OMITTED", "A1_FREQ",
                              "TEST", "OBS_CT", "BETA", "SE", "T_STAT", "P", "ERRCODE"),
             c("CHR", "BP", "SNP", "REF", "ALT", "PROVISIONAL_REF", "A1", "OMITTED", "A1_FREQ",
               "TEST", "N", "Beta", "SE", "T", "P", "ERRCODE"), skip_absent=TRUE)

    # Extract the correct interaction term for this environment
    interaction_term <- paste0("ADDx", env)  # Ensure the interaction term is correctly formatted
    plink_results <- plink_results[TEST == interaction_term]

    # Convert CHR to numeric, remove non-autosomal chromosomes
    plink_results[, CHR := as.integer(as.character(CHR))]
    plink_results <- plink_results[CHR %in% 1:22]

    # Ensure P is numeric
    plink_results[, P := as.numeric(P)]

    # Sort results
    setorder(plink_results, CHR, BP)

    # Extract filename components for saving plots
    file_name <- basename(file)
    plot_name <- gsub("\\.glm\\.linear$", "", file_name)  # Remove extension

    # Generate QQ Plot
    qqplot_file <- file.path(output_dir, paste0("qqplot_with_scaled_covar_full_dataset_", plot_name, ".png"))
    generate_qq_plot(plink_results$P, paste("QQ Plot Full Dataset:", env), qqplot_file)
  }
}

print("All QQ plots generated successfully!")
