
#clean enviroment
rm(list = ls())

#library
library(dplyr)

#set working directory
setwd("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data")

#load dataset
#load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/data2_20240930.RData")
#dataset<- dtt

# Read in the list of IDs from the file
#id_list <- read.table("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/merged_ids.rel.id", header = FALSE)$V1

# Subset the dataset to include only individuals in the ID list
#dataset <- dataset %>% filter(ID %in% id_list)  # Replace 'ID' with the actual ID column name in your dataset

#dataset$sleep_dev <- (dataset$sleep - mean(dataset$sleep))^2

#Evars <-	c("Townsend", "act0_d", "TVtime", "sleep_d", "smoking_now", "veg_cook", "fish_oily", "fish_lean", "meat_proc", "poultry", "beef", "lamb", "pork", "cheese", "salt", "tea", "alc1", "waist", "getup", "coffee", "smoked_past", "BFP", "sleep_dev")
#length(Evars)
#Emat <- as.matrix(dataset[Evars]); dim(Emat)
#rownames(Emat) <- dataset$ID
#Emat <- scale(Emat)

load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/Emat_20241106.RData")
E <- tcrossprod(Emat)
E <- E/mean(diag(E))

# Compute eigenvalues and eigenvectors
E_eigen <- eigen(E)

# Assuming row names are already set as individual IDs from 'G'
rownames(E_eigen$vectors) <- rownames(Emat)

# Save the eigenvalues and eigenvectors to an RDS file
#saveRDS(E_eigen, file = opt$output)

#save(dataset, file="data3_20241106.RData")
#save(Emat, file="Emat_20241106.RData")
#save(E, file="E_20241106.RData")
