#clean enviroment
rm(list = ls())

#library
library(dplyr)

#set working directory
setwd("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data")

#load dataset
load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/data2_20240930.RData")
dataset<- dtt

# Read in the list of IDs from the file
id_list <- read.table("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/merged_ids.rel.id", header = FALSE)$V1

# Subset the dataset to include only individuals in the ID list
dataset <- dataset %>% filter(ID %in% id_list)  # Replace 'ID' with the actual ID column name in your dataset

# Scale DP0a, PP0a, and SP0a to mean 100 and sd 10
dataset$DP0s <- scale(dataset$DP0a) * 10 + 100
dataset$PP0s <- scale(dataset$PP0a) * 10 + 100
dataset$SP0s <- scale(dataset$SP0a) * 10 + 100

# Scale AOP (age) to mean 0 and sd 1
dataset$AOPs <- scale(dataset$AOP)
dataset$AOPss<- dataset$AOPs^2

# Save the scaled dataset
save(dataset, file = "scaled_dataset_20250106.Rdata")

#set working directory
setwd("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr")

#Read in the RDS file containing the Eigen data
eigen_path_pcrelate <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/pca_for_pcrelate.rds"
eigen_results_pcrelate <- readRDS(eigen_path_pcrelate)

#Read in the RDS file containing the Eigen data
eigen_path_plink <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/pca_for_plink.rds"
eigen_results_plink <- readRDS(eigen_path_plink)

#Extract the top 10 PCs and the corresponding IDs
ids_pcrelate <- rownames(eigen_results_pcrelate$vectors)
pcs_pcrelate <- eigen_results_pcrelate$vectors[, 1:10] # Extracting the first 10 principal components

#Scale each PC to have a mean of 0 and standard deviation of 1
scaled_pcs_pcrelate <- scale(pcs_pcrelate)

#Extract the top 10 PCs and the corresponding IDs
ids_plink <- rownames(eigen_results_plink$vectors)
pcs_plink <- eigen_results_plink$vectors[, 1:10] # Extracting the first 10 principal components

#Scale each PC to have a mean of 0 and standard deviation of 1
scaled_pcs_plink <- scale(pcs_plink)

#Create a new dataset with these transformed PCs
pcs_scaled_pcrelate <- data.frame(ID = ids_pcrelate, scaled_pcs_pcrelate)

#Save the new dataset to an RDS file
saveRDS(pcs_scaled_pcrelate, file = "scaled_pcs_pcrelate.rds")

#Create a new dataset with these transformed PCs
pcs_scaled_plink <- data.frame(ID = ids_plink, scaled_pcs_plink)

#Save the new dataset to an RDS file
saveRDS(pcs_scaled_plink, file = "scaled_pcs_plink.rds")


