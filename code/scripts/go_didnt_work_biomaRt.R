rm(list=ls()); gc()
set.seed(1123)
library(data.table)
library(biomaRt)

# Use GRCh37 (hg19) for biomaRt
ensembl <- useMart(biomart = "ENSEMBL_MART_SNP", 
                   host = "https://grch37.ensembl.org", 
                   path = "/biomart/martservice", 
                   dataset = "hsapiens_snp")

# Directory containing filtered top-10 interaction term files
results_dir <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/gwas_snp_env/"
interaction_files <- list.files(results_dir, pattern="interaction_terms_top10.*\\.csv$", full.names=TRUE)

# Annotate SNPs from each file individually
for (file in interaction_files) {
  cat("Annotating file:", file, "\n")
  
  # Read file
  data <- fread(file)
  
  # Extract unique SNP IDs
  snps <- unique(data$SNP)
  
  # Run annotation using biomaRt
  annotations <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_start", 
                                      "ensembl_gene_stable_id", "associated_gene", 
                                      "consequence_type_tv"),
                       filters = "snp_filter",
                       values = snps,
                       mart = ensembl)
  
  # Generate annotation output filename
  output_filename <- gsub("interaction_terms_top10", "annotated_snps", basename(file))
  output_filepath <- file.path(results_dir, output_filename)
  
  # Save annotations
  fwrite(annotations, output_filepath)
  cat("Saved annotations to:", output_filepath, "\n")
}

cat("All files have been individually annotated and saved.\n")
