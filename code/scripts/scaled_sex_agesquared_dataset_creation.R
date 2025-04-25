##script to add scaled Sex_SI and scale(AOPss) to scaled_dataset used in analysis ##

rm(list=ls())
gc()
set.seed(1123)

#library
library(dplyr)

#set working directory
setwd("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data")

#load dataset
load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/scaled_dataset_20250106.Rdata")

dataset$AOPsss <- scale(dataset$AOPss)
dataset$Sex_SIs <- scale(dataset$Sex_SI)
length(dataset)

save(dataset, file="scaled_dataset_20250425.RData")

