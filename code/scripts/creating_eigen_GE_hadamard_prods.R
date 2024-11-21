rm(list=ls()); gc()

library(matrixcalc)
library(genio)
library(data.table)
library(Matrix)

#load GE
load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/GE_hadamard_prod_plink.RData")
GE_plink <- GE

GE_plink_eigen <- eigen(GE_plink)
# Assuming row names are already set as individual IDs from 'G'
rownames(GE_plink_eigen$vectors) <- rownames(GE_plink)

#save
saveRDS(GE_plink_eigen, file = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/eigen_GE_plink.rds")

rm(list=ls()); gc()

#load GE pcrelate
load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/GE_hadamard_prod_pcrelate.RData")
GE_pcrelate <- GE

#eigen
GE_pcrelate_eigen <- eigen(GE_pcrelate)

# Assuming row names are already set as individual IDs from 'G'
rownames(GE_pcrelate_eigen$vectors) <- rownames(GE_pcrelate)

#save
saveRDS(GE_pcrelate_eigen, file = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/eigen_GE_pcrelate.rds")


