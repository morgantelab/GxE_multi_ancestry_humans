set.seed(1123)

library(optparse)
library(data.table)

# --- Parse options ---
option_list <- list(
  make_option("--input_rds", type = "character"),
  make_option("--asian", type = "character"),
  make_option("--white", type = "character"),
  make_option("--mixed", type = "character"),
  make_option("--black", type = "character"),
  make_option("--chinese", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))

# --- Load dataset ---
dt <- readRDS(opt$input_rds)

# --- Save ancestry-specific IDs ---
ancestry_map <- list(
  ethn1_asian = opt$asian,
  ethn1_white = opt$white,
  ethn1_mixed = opt$mixed,
  ethn1_black = opt$black,
  ethn1_chinese = opt$chinese
)

for (colname in names(ancestry_map)) {
  outfile <- ancestry_map[[colname]]
  id_subset <- dt[dt[[colname]] == 1, "ID"]
  if (length(id_subset) > 0) {
    fwrite(data.table(ID1 = id_subset, ID2 = id_subset), file = outfile, sep = "\t", col.names = FALSE)
  } else {
    warning(paste("No individuals found for", colname, "- output file will be empty."))
    fwrite(data.table(ID1 = character(), ID2 = character()), file = outfile, sep = "\t", col.names = FALSE)
  }
}

cat("Saved ancestry ID files to specified paths.\n")
