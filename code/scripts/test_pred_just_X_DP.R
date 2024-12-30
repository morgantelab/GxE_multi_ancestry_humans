### Running Predictions for model with just X ###
rm(list=ls()); gc()

library(data.table)
library(Matrix)
library(BGLR)
library(readr)
library(genio)
library(dplyr)

### Please refer to https://github.com/gdlc/BGLR-R for a more valuable description ###

setwd("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model")

### load dataset ###
load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/scaled_dataset_20241025.Rdata")
print("dataset loaded")

### load fold ###
fold_1 <- readRDS("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/Fold_1.rds")
fold_1_ids <- fold_1$ID

### Extract the phenotype vectors ###
y <- dataset$DP0s
rownames(y) <- dataset$ID

print("y_PP_pcrelate created")

# Set the IDs in fold_1 to NA in y
if (!is.null(rownames(y))) {
  y[rownames(y) %in% fold_1_ids] <- NA
} else {
  stop("Row names for y are missing. Ensure y has rownames corresponding to IDs.")
}

### X is the incidence matrix for the 'fixed' covariates (no penalisation, no shrinkage). here age and sex ###
### Extract the covariate matrix ###
X <- dataset[, c("AOPs", "AOPss", "Sex_SI")]
X$AOPs <- as.vector(X$AOPs)
X$AOPss <- as.vector(X$AOPss)
rownames(X) <- dataset$ID

print("initial X created")

### Model ###

iter <- 90000
burnin <- 40000
thin <- 50
verb <- T
nrow_varabs <- (iter-burnin)/thin
scratch <- "/scratch3/kgoda/ukbiobank_files/tmp/snakemake_runs" # Use this!

### Define model ###

ETA <- list(X=list(X=X, model="FIXED", saveEffects=TRUE)
)

### Ensure ETA$X$X is numeric ###
if (!is.numeric(ETA$X$X)) {
  ETA$X$X <- as.matrix(ETA$X$X)
  ETA$X$X <- apply(ETA$X$X, 2, as.numeric)
}

print("model ETA created")

model <- BGLR(y=y, ETA=ETA, nIter=iter, burnIn=burnin, thin=thin, verbose=verb, saveAt=paste(scratch, '/pred_DP_X1_run_', sep=''))

# Extract predicted values
yHat <- model$yHat # Predicted values for the entire dataset

# Extract test set predictions
test_ids <- fold_1_ids # IDs in fold_1
if (!is.null(rownames(y))) {
  yHat_test <- yHat[rownames(y) %in% test_ids]
} else {
  stop("Row names for y are missing. Ensure y has rownames corresponding to IDs.")
}

# Combine predictions and observed values into a data frame
preds <- data.frame(ID=rownames(y), Observed=model$y, Predicted=model$yHat)

# Save predictions
results_dir <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model"

trait <- "DP0s"
analysis <- "Fold1"
save_file <- paste(results_dir, paste0('PREDS_', trait, '_BGLR_', analysis, '.RData'), sep='/')
save(preds, file=save_file)

## checking file ##
#load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/PREDS_DP0s_BGLR_Fold1.RData")
