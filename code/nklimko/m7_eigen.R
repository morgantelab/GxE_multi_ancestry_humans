.libPaths('/data2/morgante_lab/kgoda/software/R/x86_64-pc-linux-gnu-library/4.2')

# library(matrixcalc)
# library(genio)
# library(data.table)
# library(Matrix)

set.seed(123)



setwd('/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry')

m7 <- readRDS('data/nklimko/grm10.Rds')

fin7 <- eigen(m7)

saveRDS(fin7, file='data/nklimko/model7/final_eigen.rds')
