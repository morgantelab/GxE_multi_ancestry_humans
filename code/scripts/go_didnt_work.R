rm(list=ls()); gc()
set.seed(1123)
library(data.table)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(VariantAnnotation)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

# Setup genome annotation databases
snps_db <- SNPlocs.Hsapiens.dbSNP144.GRCh37
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# Directory containing filtered top-10 interaction term files
results_dir <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/gwas_snp_env/"
interaction_files <- list.files(results_dir, pattern="interaction_terms_top10.*\\.csv$", full.names=TRUE)

for (file in interaction_files) {
  cat("Annotating file:", file, "\n")
  data <- fread(file)
  snps <- unique(data$SNP)
  
  # Get genomic coordinates from rsIDs
  snp_locs <- snpsById(snps_db, snps, ifnotfound="drop")
  
  # Convert to GRanges
  snp_gr <- GRanges(seqnames = seqnames(snp_locs),
                    ranges = ranges(snp_locs),
                    rsid = snp_locs$RefSNP_id)
  
  # Annotate variants
  hits <- locateVariants(snp_gr, txdb, AllVariants())
  
  # Add gene symbols
  gene_symbols <- select(org.Hs.eg.db, keys=hits$GENEID, columns="SYMBOL", keytype="ENTREZID")
  hits_df <- as.data.table(hits)
  annotated_results <- merge(hits_df, gene_symbols, by.x="GENEID", by.y="ENTREZID", all.x=TRUE)
  
  # Clean and organize output
  annotated_results <- annotated_results[, .(rsid, seqnames, start, SYMBOL, LOCATION)]
  setnames(annotated_results, c("rsID", "Chromosome", "Position", "GeneSymbol", "VariantLocation"))
  
  # Save annotations
  output_filename <- gsub("interaction_terms_top10", "annotated_snps", basename(file))
  output_filepath <- file.path(results_dir, output_filename)
  fwrite(annotated_results, output_filepath)
  
  cat("Annotations saved to:", output_filepath, "\n")
}

cat("Annotation of all files completed successfully.\n")
