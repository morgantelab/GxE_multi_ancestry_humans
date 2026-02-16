
###Collecting results###

# mkdir -p {dp,sp,pp}/{plink,pcrelate}/model{7..8}

# Sampled regression effect
B1 <- read.table(paste0(scratch, 'ETA_X1_b.dat'), header=TRUE)
B2 <- read.table(paste0(scratch, 'ETA_X2_b.dat'), header=TRUE)
B3 <- readBinMat(paste0(scratch, 'ETA_G_b.bin'))
B4 <- readBinMat(paste0(scratch, 'ETA_E_b.bin'))
B5 <- readBinMat(paste0(scratch, 'ETA_GxE_b.bin'))
B6 <- readBinMat(paste0(scratch, 'ETA_m7_b.bin'))

# Calculate variance components
varabs <- matrix(NA, nrow_varabs, 6);
colnames(varabs) <- c("V_X1", "V_X2", "V_G", "V_E", "V_GxE", 'V_m7')

# Fill variance components
varabs[, 1] <- matrixStats::colVars(ETA$X1$X %*% t(B1))[-c(1:(burnin/thin))]
varabs[, 2] <- matrixStats::colVars(ETA$X2$X %*% t(B2))[-c(1:(burnin/thin))]
varabs[, 3] <- matrixStats::colVars(tcrossprod(ETA$G$X, B3))
varabs[, 4] <- matrixStats::colVars(tcrossprod(ETA$E$X, B4))
varabs[, 5] <- matrixStats::colVars(tcrossprod(ETA$GxE$X, B5))
varabs[, 6] <- matrixStats::colVars(tcrossprod(ETA$m7$X, B6))

# Save variance components
write.csv(varabs, file=varpath, row.names=TRUE)

