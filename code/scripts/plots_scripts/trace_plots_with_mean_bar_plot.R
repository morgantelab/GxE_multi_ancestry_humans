rm(list=ls()); gc()
set.seed(1123)

# Load necessary libraries
library(ggplot2)
library(reshape2)
library(data.table)
library(dplyr)
library(cowplot)

# Set working directory
setwd("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/plots")

# Define parameters
traits <- c("SP", "DP", "PP")
models <- c("pcrelate_GxE_pcrelate_pcs_plink_unscaled_demographics", "pcrelate_GxE_pcrelate_pcs_plink_scaled_demographics")

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

# Function to generate and save trace plots + mean bar plot
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
  
  # Melt for plotting
  varabs_melted <- melt(varabs, id.vars = "Iteration", variable.name = "Component", value.name = "Variance")
  varabs_melted2 <- varabs_melted %>% filter(Component %in% c("V_X1", "V_X2", "V_E", "V_G", "V_GxE"))
  varabs_melted2$Component <- factor(varabs_melted2$Component, levels = c("V_X1", "V_X2", "V_E", "V_G", "V_GxE"))
  
  # Mean data table for bar plot
  mean_data2 <- data.table(Component = names(varabs_means), Value = varabs_means)
  mean_data2 <- mean_data2 %>% filter(Component %in% c("V_X1", "V_X2", "V_E", "V_G", "V_GxE"))
  mean_data2$Component <- factor(mean_data2$Component, levels = c("V_X1", "V_X2", "V_E", "V_G", "V_GxE"))
  
  # Trace plot
  trace_plot <- ggplot(varabs_melted2, aes(x = Iteration, y = Variance, color = Component)) +
    geom_point(size = 0.5) +
    facet_wrap(~ Component, scales = "free_y", ncol = 2) +
    geom_hline(data = mean_data2, aes(yintercept = Value), color = "black", linetype = "dashed") +
    ggtitle(paste("Trace Plots for Variance Components -", identifier)) +
    xlab("Iteration") +
    ylab("Variance") +
    theme_minimal()
  
  # Bar plot with text labels
  bar_plot <- ggplot(mean_data2, aes(x = Component, y = Value, fill = Component)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = round(Value, 2)), vjust = -0.5, size = 4) +
    theme_minimal() +
    ylab("Mean Variance") +
    xlab("")
  
  # Combine and save
  combined_plot <- plot_grid(trace_plot, bar_plot, ncol = 1, rel_heights = c(2.5, 1.2))
  ggsave(paste0("trace_and_mean_barplot_", identifier, ".png"),
         plot = combined_plot, width = 15, height = 15)
  
}

# Run across all files
for (name in names(varabs_file_paths)) {
  generate_and_save_plots(varabs_file_paths[[name]], name)
}
