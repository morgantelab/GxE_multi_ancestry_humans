rm(list=ls()); gc()

setwd("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model")

### load dataset ###
load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/scaled_dataset_20241025.Rdata")

# Load necessary library
library(dplyr)

# Ensure ethn1_consolidated is a factor
dataset$ethn1_consolidated <- as.factor(dataset$ethn1_consolidated)

# List of dependent variables
dependent_vars <- c("DP0s", "SP0s", "PP0s")

# Initialize lists to store results
lm_models <- list()
lm_summaries <- list()

# Iterate through each dependent variable
for (y in dependent_vars) {
  cat("\nFitting linear model for:", y, "\n")
  
  # Build the formula dynamically
  formula <- as.formula(paste(y, "~ ethn1_consolidated"))
  
  # Fit the linear model
  fit <- lm(formula, data = dataset)
  
  # Save the model and summary
  lm_models[[y]] <- fit
  lm_summaries[[y]] <- summary(fit)
  
  # Print the summary
  cat("\nSummary for:", y, "\n")
  print(lm_summaries[[y]])
}

# Save models and summaries if needed
saveRDS(lm_models, "lm_models.rds")
saveRDS(lm_summaries, "lm_summaries.rds")

# Optionally, export summary coefficients to CSV for each model
for (y in dependent_vars) {
  coef_df <- as.data.frame(coef(lm_summaries[[y]]))
  write.csv(coef_df, paste0("lm_coefficients_", y, ".csv"), row.names = TRUE)
}
