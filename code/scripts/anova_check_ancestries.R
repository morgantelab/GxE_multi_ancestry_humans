rm(list=ls()); gc()

setwd("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model")

### load dataset ###
load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/scaled_dataset_20241025.Rdata")

# Load necessary library
library(dplyr)

# Ensure `ethn1_consolidated` is a factor
dataset$ethn1_consolidated <- as.factor(dataset$ethn1_consolidated)

# List of dependent variables
dependent_vars <- c("DP0s", "SP0s", "PP0s")

# Iterate through each dependent variable and perform ANOVA
results <- list()

for (y in dependent_vars) {
  cat("\nPerforming ANOVA for:", y, "\n")
  
  # Build formula dynamically
  formula <- as.formula(paste(y, "~ ethn1_consolidated"))
  
  # Perform ANOVA
  anova_model <- aov(formula, data = dataset)
  
  # Save the summary
  anova_summary <- summary(anova_model)
  results[[y]] <- anova_summary
  
  # Print the summary
  print(anova_summary)
  
  for (y in dependent_vars) {
    write.csv(as.data.frame(results[[y]][[1]]), paste0("anova_summary_", y, ".csv"))
  }
}

