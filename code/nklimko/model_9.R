set.seed(1123)

# library(argparse, lib.loc = '/data2/morgante_lab/nklimko/software/R/x86_64-pc-linux-gnu-library/4.2')

# Load required libraries
.libPaths('/data2/morgante_lab/kgoda/software/R/x86_64-pc-linux-gnu-library/4.2')

library(data.table)
library(Matrix)
library(BGLR)
library(readr)
library(genio)
library(dplyr)
library(optparse)


mainCall <- function(trait,
                     modelpath,
                     datapath,
                     eigenpath,
                     gxe_eigenpath,
                     pcpath,
                     envpath,
                     outpath,
                     varpath,
                     scratch){


  # modelpath <- 'data/nklimko/model8/eigen_irm.Rds'
  # datapath <- 'data/data3_20250514.rds'
  # eigenpath <- 'data/filtered_chr/pca_for_pcrelate.rds'
  # envpath <- 'data/E_eigen.rds'
  # pcpath <- 'data/scaled_pcs_plink.rds'
  # gxe_eigenpath <- 'data/eigen_GxE_pcrelate.rds'
  #
  # trait <- 'DP'
  stepper <- 1
  print(paste0('hark', stepper))

  # Set the working directory dynamically
  setwd('/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry')

  # Assuming the file name follows the pattern: "varabs_<type>_<grm_source>_GE.csv"
  type <- trait  # Extract the type (e.g., sp, dp, pp)

  # Load dataset
  dataset <- readRDS(datapath)
  stepper <- stepper + 1
  print(paste0('hark', stepper))
  # Load eigen results
  eigen_results <- readRDS(eigenpath)
  eigenvectors <- eigen_results$vectors
  eigenvalues <- eigen_results$values

  # Filter and scale eigenvectors by positive eigenvalues
  positive_indices <- which(eigenvalues > 0)
  filtered_eigenvectors <- eigenvectors[, positive_indices]
  filtered_eigenvalues <- eigenvalues[positive_indices]
  for(i in 1:ncol(filtered_eigenvectors)) {
    filtered_eigenvectors[, i] <- filtered_eigenvectors[, i] * sqrt(filtered_eigenvalues[i])
  }
  W <- filtered_eigenvectors
  stepper <- stepper + 1
  print(paste0('harkqwe', stepper))
  individual_ids <- rownames(W)
  matched_dataset <- dataset[match(individual_ids, dataset$ID), ]

  # Prepare phenotype vector based on type
  y <- matched_dataset[[paste0(type, "0s")]]
  stepper <- stepper + 1
  print(paste0('harkasd', stepper))
  # Prepare covariate matrix
  X <- matched_dataset[, c("AOPs", "AOPsss", "Sex_SIs")]
  X$AOPs <- as.vector(X$AOPs)
  X$AOPsss <- as.vector(X$AOPsss)
  X$Sex_SIs <- as.vector(X$Sex_SIs)
  rownames(X) <- rownames(W)

  # Load scaled PCs and match to individual IDs

  # pcs_scaled <- readRDS(pcpath)
  # matched_pcs <- pcs_scaled[match(individual_ids, pcs_scaled$ID), 2:11]
  # P <- matched_pcs

  # model8 addition
  all_eigrm <- readRDS(modelpath)
  eigenvectors <- all_eigrm$vectors
  eigenvalues <- all_eigrm$values

  positive_indices <- 1:10
  filtered_eigenvectors <- eigenvectors[, positive_indices]
  filtered_eigenvalues <- eigenvalues[positive_indices]
  for(i in 1:ncol(filtered_eigenvectors)) {
    filtered_eigenvectors[, i] <- filtered_eigenvectors[, i] * sqrt(filtered_eigenvalues[i])
  }

  add_x2 <- filtered_eigenvectors



  # Load eigen of E
  E_eigen <- readRDS(envpath)
  E_eigenvectors <- E_eigen$vectors
  E_eigenvalues <- E_eigen$values

  # Filter and scale eigenvectors by positive eigenvalues
  E_positive_indices <- which(E_eigenvalues > 0)
  E_filtered_eigenvectors <- E_eigenvectors[, E_positive_indices]
  E_filtered_eigenvalues <- E_eigenvalues[E_positive_indices]
  for(i in 1:ncol(E_filtered_eigenvectors)) {
    E_filtered_eigenvectors[, i] <- E_filtered_eigenvectors[, i] * sqrt(E_filtered_eigenvalues[i])
  }
  E <- E_filtered_eigenvectors
  add_X2_E <- E[,1:3]

  pcs_scaled_raw <- readRDS(pcpath)
  pcs_scaled <- cbind(pcs_scaled_raw, add_X2_E, add_x2)
  names(pcs_scaled) <- c('ID', paste0('PC',1:10), paste0('E',1:3), paste0('eig',1:10))

  matched_pcs <- pcs_scaled[match(individual_ids, pcs_scaled$ID), 2:24]
  X2mod <- matched_pcs


  # stepper <- stepper + 1
  # print(paste0('harkdfg', stepper))
  # # Model 7 addition ----
  # print('MINGUS 7')
  # # Load eigen of E
  # E_eigen <- readRDS(modelpath)
  # E_eigenvectors <- E_eigen$vectors
  # E_eigenvalues <- E_eigen$values
  #
  # # Filter and scale eigenvectors by positive eigenvalues
  # E_positive_indices <- which(E_eigenvalues > 0)
  #
  # print('vvvv m7 positive indices vvvv ')
  # print(length(E_positive_indices))
  #
  # E_filtered_eigenvectors <- E_eigenvectors[, E_positive_indices]
  # E_filtered_eigenvalues <- E_eigenvalues[E_positive_indices]
  # for(i in 1:ncol(E_filtered_eigenvectors)) {
  #   E_filtered_eigenvectors[, i] <- E_filtered_eigenvectors[, i] * sqrt(E_filtered_eigenvalues[i])
  # }
  #
  # m7 <- E_filtered_eigenvectors
  #
  # print('m7 done')
  # Load eigen of GE
  GE_eigen <- readRDS(gxe_eigenpath)
  GE_eigenvectors <- GE_eigen$vectors
  GE_eigenvalues <- GE_eigen$values

  # Filter and scale eigenvectors by positive eigenvalues
  GE_positive_indices <- which(GE_eigenvalues > 0)
  GE_filtered_eigenvectors <- GE_eigenvectors[, GE_positive_indices]
  GE_filtered_eigenvalues <- GE_eigenvalues[GE_positive_indices]
  for(i in 1:ncol(GE_filtered_eigenvectors)) {
    GE_filtered_eigenvectors[, i] <- GE_filtered_eigenvectors[, i] * sqrt(GE_filtered_eigenvalues[i])
  }
  GE <- GE_filtered_eigenvectors

  # Model setup
  iter <- 90000
  burnin <- 40000
  thin <- 50
  verb <- T
  nrow_varabs <- (iter-burnin)/thin
  print('etamake')
  ETA <- list(X1=list(X=X, model="FIXED", saveEffects=TRUE),
              X2=list(X=X2mod, model="FIXED", saveEffects=TRUE),
              G=list(X=W, model="BRR", saveEffects=TRUE),
              E=list(X=E, model="BRR", saveEffects=TRUE),
              GxE=list(X=GE, model="BRR", saveEffects=TRUE)
  )

  if (!is.numeric(ETA$X1$X)) {
    ETA$X1$X <- as.matrix(ETA$X1$X)
    ETA$X1$X <- apply(ETA$X1$X, 2, as.numeric)
  }
  if (!is.numeric(ETA$X2$X)) {
    ETA$X2$X <- as.matrix(ETA$X2$X)
    ETA$X2$X <- apply(ETA$X2$X, 2, as.numeric)
  }
  print("Model ETA created")

  # Run BGLR model
  model <- BGLR(y=y, ETA=ETA, nIter=iter, burnIn=burnin, thin=thin, verbose=verb, saveAt=scratch)


  print(model)



  saveRDS(model, file=outpath)


  ###Collecting results###

  # mkdir -p {dp,sp,pp}/{plink,pcrelate}/model{7..8}


  print('bin read all')
  # Sampled regression effect
  B1 <- read.table(paste0(scratch, 'ETA_X1_b.dat'), header=TRUE)
  B2 <- read.table(paste0(scratch, 'ETA_X2_b.dat'), header=TRUE)
  B3 <- readBinMat(paste0(scratch, 'ETA_G_b.bin'))
  B4 <- readBinMat(paste0(scratch, 'ETA_E_b.bin'))
  B5 <- readBinMat(paste0(scratch, 'ETA_GxE_b.bin'))
  # B6 <- readBinMat(paste0(scratch, 'ETA_m8_b.bin'))

  # Calculate variance components
  varabs <- matrix(NA, nrow_varabs, 5);
  colnames(varabs) <- c("V_X1", "V_X2", "V_G", "V_E", "V_GxE")

  print('fill variance')

  # Fill variance components
  varabs[, 1] <- matrixStats::colVars(ETA$X1$X %*% t(B1))[-c(1:(burnin/thin))]
  varabs[, 2] <- matrixStats::colVars(ETA$X2$X %*% t(B2))[-c(1:(burnin/thin))]
  varabs[, 3] <- matrixStats::colVars(tcrossprod(ETA$G$X, B3))
  varabs[, 4] <- matrixStats::colVars(tcrossprod(ETA$E$X, B4))
  varabs[, 5] <- matrixStats::colVars(tcrossprod(ETA$GxE$X, B5))
  # varabs[, 6] <- matrixStats::colVars(tcrossprod(ETA$m8$X, B6))

  # Save variance components
  write.csv(varabs, file=varpath, row.names=TRUE)


}

stepper <- 1
print(paste0('hark', stepper))
print('pastmain')
#
# #args----
# parser <- ArgumentParser(description= 'snakemake transfer')
# stepper <- stepper + 1
# print(paste0('hark', stepper))
# parser$add_argument('--trait')
# parser$add_argument('--modelpath')
# parser$add_argument('--datapath')
# parser$add_argument('--eigenpath')
# parser$add_argument('--gxe_eigenpath')
# parser$add_argument('--pcpath')
# parser$add_argument('--envpath')
# parser$add_argument('--outpath')
# parser$add_argument('--varpath')
# parser$add_argument('--scratch')
# stepper <- stepper + 1
# print(paste0('hark', stepper))
# snake <- parser$parse_args()
# stepper <- stepper + 1
# print(paste0('hark', stepper))
# print(str(snake))
# stepper <- stepper + 1
# print(paste0('hark', stepper))
# print('about to call')
#
#
#

option_list <- list(make_option(c('--trait', type='character')),
                    make_option(c('--modelpath', type='character')),
                    make_option(c('--datapath', type='character')),
                    make_option(c('--eigenpath', type='character')),
                    make_option(c('--gxe_eigenpath', type='character')),
                    make_option(c('--pcpath', type='character')),
                    make_option(c('--envpath', type='character')),
                    make_option(c('--outpath', type='character')),
                    make_option(c('--varpath', type='character')),
                    make_option(c('--scratch', type='character')))



opt_parser <- OptionParser(option_list = option_list)
snake <- parse_args(opt_parser)


#call----
mainCall(snake$trait,
         snake$modelpath,
         snake$datapath,
         snake$eigenpath,
         snake$gxe_eigenpath,
         snake$pcpath,
         snake$envpath,
         snake$outpath,
         snake$varpath,
         snake$scratch)

stepper <- stepper + 1
print(paste0('harkonen', stepper))


