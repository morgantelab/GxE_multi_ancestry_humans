rm(list=ls()); gc()

# Set seed for reproducibility
set.seed(1123)

library(data.table)

setwd("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/gwas_snp_env/")
data <- fread("covar_pheno.txt")

# Scale the 'sex' and 'age2' columns
data[, `:=`(sex = scale(sex), age2 = scale(age2))]

# Save the modified dataset
fwrite(data, "covar_pheno_sex_age2_scaled.txt", sep="\t", quote=FALSE)

