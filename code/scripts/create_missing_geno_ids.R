set.seed(1123)

suppressWarnings(suppressMessages({
  library(optparse)
  library(data.table)
}))

# --- Parse command-line arguments ---
option_list <- list(
  make_option("--rds_file", type = "character", help = "Path to .rds file with IDs"),
  make_option("--fam_files", type = "character", action = "store", help = "Paths to .fam files", nargs = Inf),
  make_option("--out_file", type = "character", help = "Path to output file for missing IDs")
)
opt <- parse_args(OptionParser(option_list = option_list))

# --- Load dataset ---
dt <- readRDS(opt$rds_file)
all_ids <- unique(as.character(dt$ID))

# --- Collect present IDs from all .fam files ---
present_ids <- c()
for (fam in opt$fam_files) {
  if (!file.exists(fam)) {
    warning(paste("Missing .fam file:", fam))
    next
  }
  fam_data <- fread(fam, header = FALSE)
  present_ids <- c(present_ids, fam_data$V2)
}
present_ids <- unique(as.character(present_ids))

# --- Compute and write missing IDs ---
missing_ids <- sort(setdiff(all_ids, present_ids))
writeLines(missing_ids, con = opt$out_file)
cat(length(missing_ids), "individuals found with missing genotype data.\n")
