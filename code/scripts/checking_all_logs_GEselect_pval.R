rm(list=ls()); gc()
set.seed(1123)

# Set the path to your logs directory
log_dir <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/output/logs_slurm"

# List all .err files in the directory
all_err_files <- list.files(log_dir, pattern = "\\.err$", full.names = TRUE)

# Define selection criteria
traits <- c("SP", "DP", "PP")
folds <- 1:5
models <- c("X1_X2_G_E_GE_GEselect")
pvals <- c("1e-05", "1e-06", "0.001")
ancestries <- c("asian", "mixed", "white", "black", "chinese")  # allowed ancestries

# Build regex patterns:
# For fold-based logs, e.g., 
# "pred_model_fold_X1_X2_G_E_GE_GEselect_bp=DP,foldnum=3,pval=0.001..."
fold_pattern <- paste0("pred_model_fold_", models, "_bp=(",
                       paste(traits, collapse = "|"), "),foldnum=([1-5]),pval=(",
                       paste(pvals, collapse = "|"), ")")

# For ethnicity-based logs, e.g., 
# "pred_model_ethn_X1_X2_G_E_GE_GEselect_ancestry=white,bp=DP,pval=1e-06..."
ethn_pattern <- paste0("pred_model_ethn_", models, "_ancestry=(",
                       paste(ancestries, collapse = "|"), "),bp=(",
                       paste(traits, collapse = "|"), "),pval=(",
                       paste(pvals, collapse = "|"), ")")

# Filter files based on these patterns
fold_files <- grep(fold_pattern, all_err_files, value = TRUE)
ethn_files <- grep(ethn_pattern, all_err_files, value = TRUE)

# Combine the selected files
selected_err_files <- c(fold_files, ethn_files)

# Optional: print selected file names for verification
cat("Selected .err files:\n")
print(selected_err_files)

# For each selected file, check if it contains the string "Iter=90000"
iteration_check <- sapply(selected_err_files, function(file) {
  lines <- readLines(file, warn = FALSE)
  any(grepl("Iter=90000", lines))
})

# Create a summary data frame
results_df <- data.frame(
  file = selected_err_files,
  reached_90000 = iteration_check,
  stringsAsFactors = FALSE
)

# Display the results
print(results_df)
