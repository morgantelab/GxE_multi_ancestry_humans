rm(list=ls()); gc()
library(data.table)
library(qqman)
library(ggplot2)
library(genio)  # Read PLINK files for SNP metadata

# Define file paths
plink_results_file <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/gwas_snp_env/gxe_meat_proc_1.DP0s.glm.linear"
plink_prefix <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/merged_geno_common_snps_selected_indiv"

# Load PLINK GWAS results
plink_results <- fread(plink_results_file)

# Rename columns based on actual PLINK output
setnames(plink_results, c("#CHROM", "POS", "ID", "REF", "ALT", "PROVISIONAL_REF?", "A1", "OMITTED", "A1_FREQ",
                          "TEST", "OBS_CT", "BETA", "SE", "T_STAT", "P", "ERRCODE"),
         c("CHR", "BP", "SNP", "REF", "ALT", "PROVISIONAL_REF", "A1", "OMITTED", "A1_FREQ",
           "TEST", "N", "Beta", "SE", "T", "P", "ERRCODE"))

# Keep only the SNP × Environment interaction term
plink_results <- plink_results[TEST == "ADDxbeef"]  # Adjust `TEST` based on your interaction term

# Convert CHR to numeric, remove non-autosomal chromosomes
plink_results[, CHR := as.integer(as.character(CHR))]
plink_results <- plink_results[CHR %in% 1:22]  # Keep only chromosomes 1-22

# Load SNP metadata using `genio`
plink_data <- read_plink(plink_prefix)  # Reads .bim, .bed, .fam files
snp_metadata <- as.data.table(plink_data$bim)  # Convert bim file to data.table

# Ensure SNP metadata has correct column names
setnames(snp_metadata, c("CHR", "SNP", "CM", "BP", "A1", "A2"))

# Merge SNP metadata with GWAS results
plink_results <- merge(plink_results, snp_metadata[, c("CHR", "SNP", "BP")], by="SNP", all.x=TRUE)

# Fix CHR column naming issue after merge
if ("CHR.x" %in% colnames(plink_results)) {
  setnames(plink_results, c("CHR.x", "BP.x"), c("CHR", "BP"))
}

# Remove duplicate CHR and BP columns from merge (e.g., CHR.y, BP.y)
plink_results <- plink_results[, !grepl("\\.y$", names(plink_results)), with=FALSE]


# Remove missing values and ensure correct data types
plink_results <- plink_results[!is.na(CHR) & !is.na(BP)]
plink_results[, P := as.numeric(P)]

# Sort results for Manhattan plotting
setorder(plink_results, CHR, BP)

# **Manhattan Plot for SNP × beef Interaction**
png("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/plots/manhattan_plink_gxe_beef_1_DP.png", width=1200, height=600)
manhattan(plink_results, chr="CHR", bp="BP", p="P", snp="SNP",
          main="Manhattan Plot: PLINK SNP × beef Interaction",
          ylim=c(0, -log10(min(plink_results$P, na.rm=TRUE)) + 1),
          col=c("blue4", "orange3"), suggestiveline=-log10(1e-5), genomewideline=-log10(5e-8))
dev.off()

# **QQ Plot for SNP × Environment Interaction**
png("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/plots/qqplot_plink_gxe_beef_1_DP.png", width=800, height=600)
qq(plink_results$P, main="QQ Plot: PLINK SNP × beef Interaction")
dev.off()

print("Plots saved successfully!")
