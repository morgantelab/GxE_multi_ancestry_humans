rm(list=ls()); gc()
# Load necessary libraries
library(SNPRelate)
library(GENESIS)
library(GWASTools)
library(BiocParallel)

# Load genotype data for PCA analysis
diversity_geno <- GdsGenotypeReader(filename = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/merged_geno_common_snps_selected_indiv.gds")
diversity_genoData <- GenotypeData(diversity_geno)

# Load PC-AiR results for PC-Relate
diversity_genoData <- GenotypeBlockIterator(diversity_genoData)

#load pcair_results
pcair_results <- readRDS("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/pcair_for_pcrelate_results.rds")

# Run PC-Relate with 5 principal components
mypcrelate_5pcs <- pcrelate(
  diversity_genoData,
  pcs = pcair_results$vectors[, 1:5],  # Use first 5 PCs
  ibd.probs = FALSE,
  training.set = pcair_results$unrels,
  BPPARAM = BiocParallel::MulticoreParam(workers = 55)
)

# Save the PC-Relate results
saveRDS(mypcrelate_5pcs, file = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/test_mypcrelate_R_412.RDS")

print("Analysis completed.")

