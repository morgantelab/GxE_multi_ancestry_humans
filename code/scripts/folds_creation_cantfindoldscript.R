# Load the dataset
load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/scaled_dataset_20241025.Rdata")

# Assuming your dataset is stored in a variable named `data`
# If not, replace `data` with the actual variable name in the loaded file.

# Set seed for reproducibility
set.seed(123)

# Check the number of rows
n <- nrow(data)

# Create a Fold column with random assignment
data$Fold <- sample(rep(1:5, length.out = n))

# Check the distribution of folds
table(data$Fold)

# Split the dataset into 5 folds if needed
folds <- split(data, data$Fold)

# Save each fold separately if required
for (i in 1:5) {
  saveRDS(folds[[i]], paste0("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/Fold_", i, ".Rds"))
}

# Optional: View the structure of the folds
str(folds)
