rm(list=ls()); gc()
set.seed(1123)

# Load necessary libraries
library(ggplot2)
library(reshape2)
library(data.table)
library(dplyr)

# Set working directory
setwd("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/plots")

# Define parameters
traits <- c("SP", "DP", "PP")
#models <- c("S_A_X2", "S_A_X2_G", "S_A_X2_E", "S_A_X2_G_E", "S_A_X2_G_E_GE")
models <- c("pcrelate_GxE_pcrelate_pcs_plink_unscaled_demographics")

# Generate filenames
varabs_file_paths <- setNames(
  unlist(lapply(traits, function(trait) {
    lapply(models, function(model) {
      paste0("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_", trait, "_", model, ".csv")
    })
  })),
  unlist(lapply(traits, function(trait) {
    lapply(models, function(model) {
      paste0(trait, "_", model)
    })
  }))
)

# Function to generate and save trace plots
generate_and_save_plots <- function(varabs_file_path, identifier) {
  if (!file.exists(varabs_file_path)) {
    message("File not found: ", varabs_file_path)
    return()
  }

  # Read varabs data
  varabs <- fread(varabs_file_path)

  # Convert all columns to double, if possible
  varabs <- varabs[, lapply(.SD, as.double)]

  # Calculate the mean for each component
  varabs_means <- colMeans(varabs, na.rm = TRUE)

  # Add iteration column
  varabs$Iteration <- seq(1, nrow(varabs))

  # Melt data for plotting
  varabs_melted <- melt(varabs, id.vars = "Iteration", variable.name = "Component", value.name = "Variance")
  varabs_melted2 <- varabs_melted %>% filter(Component %in% c("V_X1" ,"V_X2", "V_E", "V_G", "V_GxE"))

  varabs_melted2$Component <- factor(varabs_melted2$Component, levels = c("V_X1", "V_X2", "V_E", "V_G", "V_GxE"))

  # Create a mean data table
  mean_data2 <- data.table(Component = names(varabs_means), Value = varabs_means)
  mean_data2 <- mean_data2 %>% filter(Component %in% c("V_X1", "V_X2", "V_E", "V_G", "V_GxE"))

  # Plot trace plots
  variances_plot <- ggplot(varabs_melted2, aes(x = Iteration, y = Variance, color = Component)) +
    geom_point(size = 0.5) +
    facet_wrap(~ Component, scales = "free_y", ncol = 2) +
    geom_hline(data = mean_data2, aes(yintercept = Value), color = "black", linetype = "dashed") +
    ggtitle(paste("Trace Plots for Variance Components -", identifier)) +
    xlab("Iteration") +
    ylab("Variance")

  # Save plot
  ggsave(paste0("trace_plots_variance_components_", identifier, ".png"), plot = variances_plot, width = 10, height = 6)
}

# Iterate over each generated file path
for (name in names(varabs_file_paths)) {
  generate_and_save_plots(varabs_file_paths[[name]], name)
}
