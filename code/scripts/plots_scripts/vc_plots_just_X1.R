rm(list=ls()); gc()

# Load necessary libraries
library(ggplot2)
library(reshape2)
library(data.table)
library(dplyr)

# Set working directory
setwd("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/plots")

# Define lists of file paths for varabs
varabs_file_paths <- list(
  DP_pcrelate = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_DP_X.csv",
  SP_pcrelate = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_SP_X.csv",
  PP_pcrelate = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_PP_X.csv"
)

# Define function to generate and save plots for each file
generate_and_save_plots <- function(varabs_file_path, identifier) {
  
  # Read the varabs data
  varabs <- fread(varabs_file_path)
  
  # Convert all columns to double, if possible
  varabs <- varabs[, lapply(.SD, as.double)]
  
  # Calculate the mean for each component
  varabs_means <- colMeans(varabs, na.rm = TRUE)
  
  # Add iteration column
  varabs$Iteration <- seq(1, nrow(varabs))
  
  # Melt the varabs data for plotting
  varabs_melted <- melt(varabs, id.vars = "Iteration", variable.name = "Component", value.name = "Variance")
  
  varabs_melted2 <- varabs_melted %>% filter(!Component == "V1")
  varabs_melted2$Component <- factor(varabs_melted2$Component, levels = c("V_X"))
  
  # Create a mean data table with the same factor levels
  mean_data2 <- data.table(Component = factor(names(varabs_means), levels = levels(varabs_melted$Component)),
                           Value = varabs_means)
  mean_data2 <- mean_data2 %>% filter(!Component == "V1")
  
  
  # Plot trace plots for varabs components
  variances_plot <- ggplot(varabs_melted2, aes(x = Iteration, y = Variance, color = Component)) +
    geom_point(size = 0.5) +
    facet_wrap(~ Component, scales = "free_y") +
    geom_hline(data = mean_data2, aes(yintercept = Value), color = "black", linetype = "dashed") +
    ggtitle(paste("Trace Plots for Variance Components -", identifier)) +
    xlab("Iteration") +
    ylab("Variance")
  
  # Save the varabs plot
  ggsave(paste0("with_mean_X_trace_plots_variance_components_", identifier, ".png"), plot = variances_plot, width = 10, height = 6)
}

# Iterate over each pair of VCEm and varabs file paths
for (name in names(varabs_file_paths)) {
  varabs_file_path <- varabs_file_paths[[name]]
  generate_and_save_plots(varabs_file_path, name)
}

