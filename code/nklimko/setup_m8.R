.libPaths('/data2/morgante_lab/kgoda/software/R/x86_64-pc-linux-gnu-library/4.2')

library(matrixcalc)
library(genio)
library(data.table)
library(Matrix)



setwd('/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry')

# Creates gcovariates for model 8 ----
if(1) {

# erm creation ----
  not_made <- FALSE
  if(not_made){
    # load env data, 24k x 23 matrix
    load('data/Emat_20250514.RData')

    # Compute E from Emat
    # {E * E^T} / mean(diag(E)
    E <- tcrossprod(Emat)
    etude <- E / mean(diag(E))

    saveRDS(etude, file='data/nklimko/erm.Rds')
  } else {
    # etude <- readRDS('data/nklimko/erm.Rds')
  }



  ## Interaction relation matrix creation----
  # not_made <- FALSE
  # if(not_made){
  #   plink_grm_data <- read_grm('data/filtered_chr/grm_by_plink')
  #   G_plink <- plink_grm_data$kinship
  #   G_plink <- as.matrix(G_plink)
  #
  #   irm <- etude * G_plink
  #   saveRDS(irm, file='data/nklimko/irm.Rds')
  # } else {
  #   irm <- readRDS('data/nklimko/irm.Rds')
  # }



  # Eigen decompositions----
  not_made <- TRUE
  if(not_made){
    irm_eigen <- eigen(irm)
    rownames(irm_eigen$vectors) <- rownames(irm)
    saveRDS(irm_eigen, file = 'data/nklimko/model8/eigen_irm.Rds')

  } else {
    irm_eigen <- readRDS('data/nklimko/model8/eigen_irm.Rds')
  }


  # # subset 10 covariates ----
  # elect <- 10
  # vec10 <- irm_eigen$vectors[,1:elect]
  # saveRDS(vec10, file='data/nklimko/model8/eigvec10.rds')

  # val10 <- irm_eigen$values[1:elect]
  # tvec10 <- t(vec10)
  # dag10 <- diag(val10)
  # g10 <- vec10 %*% dag10 %*% tvec10

  # ud <- vec10 %*% dag10
  # udut <- ud %*% tvec10

}



# ws ----
# all_eigrm <- readRDS('data/nklimko/model8/eigen_irm.Rds')
#
#
#   all_eigrm <- readRDS(modelpath)
#   eigenvectors <- all_eigrm$vectors
#   eigenvalues <- all_eigrm$values
#
#   positive_indices <- 1:10
#   filtered_eigenvectors <- eigenvectors[, positive_indices]
#   filtered_eigenvalues <- eigenvalues[positive_indices]
#   for(i in 1:ncol(filtered_eigenvectors)) {
#     filtered_eigenvectors[, i] <- filtered_eigenvectors[, i] * sqrt(filtered_eigenvalues[i])
#   }
#
#   add_x2 <- filtered_eigenvectors
#
#
#     pcs_scaled_raw <- readRDS(pcpath)
#
#     pcs_scaled <- cbind(pcs_scaled_raw, add_x2)
#
#     names(pcs_scaled) <- c('ID', paste0('PC',1:10),  paste0('eig',1:10))
#
#   matched_pcs <- pcs_scaled[match(individual_ids, pcs_scaled$ID), 2:21]
#   X2mod <- matched_pcs
#
# dim(W)
# W
# rownames(W)
#
#
# pcpath <- 'data/scaled_pcs_plink.rds'
#
#
# pcs_scaled$ID[1:5]
#
# rownames(W)[1:5]
#
# sum(pcs_scaled$ID == rownames(W))

