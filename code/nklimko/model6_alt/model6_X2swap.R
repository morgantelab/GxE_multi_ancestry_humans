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
                     datapath,
                     eigenpath,
                     gxe_eigenpath,
                     pcpath,
                     envpath,
                     outpath,
                     varpath,
                     scratch){


  stepper <- 1
  # print(paste0('hark', stepper))

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
  pcs_scaled <- readRDS(pcpath)
  matched_pcs <- pcs_scaled[match(individual_ids, pcs_scaled$ID), 2:11]
  P <- matched_pcs
  stepper <- stepper + 1
  print(paste0('harksdf', stepper))
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
  # stepper <- stepper + 1
  # print(paste0('harkdfg', stepper))

  # Model 7 addition ----
  print('MINGUS 7')
  # Load eigen of E
  if(FALSE){
    E_eigen <- readRDS(modelpath)
    E_eigenvectors <- E_eigen$vectors
    E_eigenvalues <- E_eigen$values

    # Filter and scale eigenvectors by positive eigenvalues
    E_positive_indices <- which(E_eigenvalues > 0)

    print('vvvv m7 positive indices vvvv ')
    print(length(E_positive_indices))

    E_filtered_eigenvectors <- E_eigenvectors[, E_positive_indices]
    E_filtered_eigenvalues <- E_eigenvalues[E_positive_indices]
    for(i in 1:ncol(E_filtered_eigenvectors)) {
      E_filtered_eigenvectors[, i] <- E_filtered_eigenvectors[, i] * sqrt(E_filtered_eigenvalues[i])
    }
    m7 <- E_filtered_eigenvectors
  }

  # m7 <- readRDS(modelpath)

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
              X2=list(X=P, model="FIXED", saveEffects=TRUE),
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
  # B6 <- readBinMat(paste0(scratch, 'ETA_M7_b.bin'))

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
  # varabs[, 6] <- matrixStats::colVars(tcrossprod(ETA$M7$X, B6))

  # Save variance components
  write.csv(varabs, file=varpath, row.names=TRUE)


}

stepper <- 1
# print(paste0('hark', stepper))
print('pastmain')

option_list <- list(make_option(c('--trait', type='character')),
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


