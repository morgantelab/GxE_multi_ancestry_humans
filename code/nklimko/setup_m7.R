library(matrixcalc)
library(genio)
library(data.table)
library(Matrix)
library(optparse)

# Creates g10 subset for model 7 ----
if(1) {

  # top 10 pcs
  elect <- 10

  if(0){
    pca_plink <- readRDS('data/filtered_chr/pca_for_plink.rds')
    val10 <- pca_plink$values[1:elect]
    vec10 <- pca_plink$vectors[,1:elect]
  }

  # pca pcrelate ----
  pca_pcrelate <- readRDS('data/filtered_chr/pca_for_pcrelate.rds')
  val10 <- pca_pcrelate$values[1:elect]
  vec10 <- pca_pcrelate$vectors[,1:elect]

  tvec10 <- t(vec10)
  dag10 <- diag(val10)

  # g10 <- vec10 %*% dag10 %*% tvec10

  ud <- vec10 %*% dag10
  udut <- ud %*% tvec10

  saveRDS(udut, file='data/nklimko/g10.rds')

  # ERM ----
  load('data/Emat_20250514.RData')
  # dim(Emat)

  # Compute E from Emat
  E <- tcrossprod(Emat)

  etude <- E / mean(diag(E))
  # dim(etude)

  saveRDS(etude, file='data/nklimko/erm.Rds')
  # dim(udut)

# IRM ----
  etude <- readRDS('data/nklimko/erm.Rds')
  m7 <- udut * etude

  saveRDS(m7, file='data/nklimko/m7_pcrelate_irm10.rds')

  # eigen decomp ----
  m7 <- readRDS('data/nklimko/grm10.Rds')
  fin7 <- eigen(m7)

  saveRDS(fin7, file='data/nklimko/model7/final_eigen.rds')
  # ws ----
  # pcs <- prcomp(m7, scale=TRUE)

  # saveRDS(pcs, file='data/nklimko/model7/allpcs.rds')

  # pc10 <- pcs[,1:10]

  # saveRDS(pc10, file='data/nklimko/model7/pc10.rds')

}



# eigen script, split for runtime speedup in m7_eigen.r ----
# .libPath('/data2/morgante_lab/kgoda/software/R/x86_64-pc-linux-gnu-library/4.2')
#
# library(matrixcalc)
# library(genio)
# library(data.table)
# library(Matrix)
#
# set.seed(123)
#
# setwd('/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry')
# m7 <- readRDS('data/nklimko/grm10.Rds')
# fin7 <- eigen(m7)
# saveRDS(fin7, file='data/nklimko/model7/final_eigen.rds')


# ws ----
# temp <- readRDS('data/nklimko/grm10.Rds')
# x2_pcs <- readRDS('data/scaled_pcs_plink.rds')
# dim(x2_pcs)
# m7
