rm(list=ls()); gc()
set.seed(1123)

library(optparse)
library(data.table)

# Define command-line options
option_list <- list(
  make_option(c("-i", "--input_dir"), type = "character", default = NULL,
              help = "Directory containing PLINK results", metavar = "character"),
  make_option(c("-o", "--output_dir"), type = "character", default = NULL,
              help = "Directory to save filtered results", metavar = "character"),
  make_option(c("-p", "--pval_threshold"), type = "numeric", default = 1e-5,
              help = "P-value threshold for filtering", metavar = "numeric")
)

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$input_dir) || is.null(opt$output_dir)) {
  stop("Both input and output directories must be specified.")
}

# Define known interaction terms
interaction_terms <- c("Townsend", "act0_d", "TVtime", "sleep_d", "smoking_now", "veg_cook", "fish_oily",
                       "fish_lean", "meat_proc", "poultry", "beef", "lamb", "pork", "cheese", "salt", "tea",
                       "alc1", "waist", "getup", "coffee", "smoked_past", "BFP", "sleep_dev")

# List all PLINK result files
result_files <- list.files(path = opt$input_dir, pattern = "*.glm.linear$", full.names = TRUE)

# Initialize summary table
summary_results <- data.table(File = character(), Hits = integer(), Total_Tests = integer(), Mean_P_Filtered = numeric(), Mean_P_All = numeric())

for (interaction_term in interaction_terms) {
  print(paste("🔹 Processing interaction:", interaction_term))

  selected_files <- result_files[grepl(paste0("_", interaction_term, "_"), result_files)]

  if (length(selected_files) == 0) {
    print(paste("⚠️ No files found for interaction:", interaction_term))
    next
  }

  for (plink_results_file in selected_files) {
    file_name <- basename(plink_results_file)
    print(paste("✅ Processing file:", file_name, "for interaction term:", interaction_term))

    results <- fread(plink_results_file)

    if (!("TEST" %in% colnames(results)) || !("P" %in% colnames(results))) {
      print(paste("⚠️ Skipping", file_name, "- required columns (TEST, P) not found"))
      next
    }

    results[, P := as.numeric(P)]
    results <- results[!is.na(P)]

    interaction_test <- paste0("ADDx", interaction_term)
    results_interaction <- results[TEST == interaction_test]
    results_filtered <- results_interaction[P < opt$pval_threshold]

    num_hits <- nrow(results_filtered)
    total_tests <- nrow(results_interaction)
    mean_p_filtered <- ifelse(num_hits > 0, mean(results_filtered$P), NA)
    mean_p_all <- ifelse(total_tests > 0, mean(results_interaction$P), NA)

    summary_results <- rbind(summary_results, data.table(File = file_name, Hits = num_hits, Total_Tests = total_tests, Mean_P_Filtered = mean_p_filtered, Mean_P_All = mean_p_all))

    if (num_hits == 0) {
      print(paste("⚠️ No significant hits found in", file_name))
      next
    }

    output_file <- file.path(opt$output_dir, paste0(sub(".glm.linear", "", file_name), "_filtered_", opt$pval_threshold, ".txt"))
    fwrite(results_filtered, output_file, sep = "\t", quote = FALSE)
    print(paste("💾 Saved filtered results to", output_file))
  }
}

# Save summary results
summary_file <- file.path(opt$output_dir, paste0("summary_hits_", opt$pval_threshold, ".txt"))
fwrite(summary_results, summary_file, sep = "\t", quote = FALSE)
print(paste("📊 Summary of hits saved to", summary_file))

print("✅ Processing complete.")
