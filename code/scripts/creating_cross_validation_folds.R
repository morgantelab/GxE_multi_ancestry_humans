##Creating 5-folds for the dataset##

#clear working space
rm(list=ls()); gc()

#set working directory
setwd("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model")

# Load the dataset
load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/scaled_dataset_20241025.Rdata")

# Set seed for reproducibility
set.seed(123)

# Check the number of rows
n <- nrow(dataset)

# Create a Fold column with random assignment
dataset$Fold <- sample(rep(1:5, length.out = n))

# Check the distribution of folds
table(dataset$Fold)

# Split the dataset into 5 folds if needed
folds <- split(dataset, dataset$Fold)

# Save each fold separately if required
for (i in 1:5) {
  saveRDS(folds[[i]], paste0("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/Fold_", i, ".rds"))
}

# Optional: View the structure of the folds
str(folds)
