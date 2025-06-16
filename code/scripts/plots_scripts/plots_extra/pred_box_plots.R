rm(list=ls()); gc()

# Load necessary library
library(ggplot2)

setwd("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/plots")

# Load the dataset
dataset <- read_csv("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/final pred numbers.csv")

# Display the first few rows of the dataset to understand its structure
head(dataset)

# Define the output directory
output_dir <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/plots"

# Function to create and save plots for a given metric (e.g., "R^2" or "Correlation")
create_and_save_plots <- function(metric, metric_label) {
  for (trait in unique(dataset$Trait)) {
    # Filter the data for the specific trait
    trait_data <- subset(dataset, Trait == trait)
    
    # Create the box plot
    plot <- ggplot(trait_data, aes(x = Model, y = get(metric))) +
      geom_boxplot() +
      labs(
        title = paste("Box Plot of", metric_label, "for", trait, "by Model"),
        x = "Model",
        y = metric_label
      ) +
      theme_minimal()
    
    # Save the plot with a unique file name
    file_name <- paste(output_dir, "/boxplot_", trait, "_", metric, ".png", sep = "")
    ggsave(file_name, plot = plot)
  }
}

# Create and save plots for R^2
create_and_save_plots("R^2", "R-squared")

# Create and save plots for Correlation
create_and_save_plots("Correlation", "Correlation")

# Confirm completion
print(paste("All plots for R^2 and Correlation have been saved to", output_dir))


# Function to create styled box plots with outlined data points
create_styled_boxplot <- function(data, y_variable, y_label, title, file_name) {
  plot <- ggplot(data, aes(x = Model, y = !!sym(y_variable), fill = Model)) +
    geom_boxplot(outlier.shape = NA, color = "black", width = 0.6, lwd = 0.8) +  # Thicker borders
    geom_point(aes(group = Model), shape = 1, size = 3, stroke = 1.2, color = "black", position = position_jitter(width = 0.15)) +  # Outlined points
    scale_fill_manual(values = c("#F4A582", "#92C5DE", "#D6604D", "#4393C3", "#FDAE61", "#67A9CF")) + # Custom colors
    labs(
      title = title,
      x = "Model",
      y = y_label
    ) +
    theme_classic(base_size = 14) +  # Clean theme with customizations
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),  # Rotated labels
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.position = "none",  # Remove legend
      axis.line = element_line(size = 0.8, color = "black")
    )
  
  # Save the plot
  ggsave(
    filename = paste0(output_dir, "/", file_name),
    plot = plot,
    width = 8, height = 6, dpi = 300
  )
}

# Loop through each trait and create plots for both R^2 and Correlation
traits <- unique(dataset$Trait)
for (trait in traits) {
  # Filter the data for the specific trait
  trait_data <- subset(dataset, Trait == trait)
  
  # Create R^2 plot
  create_enhanced_plot(
    data = trait_data,
    y_variable = "R^2",
    y_label = expression(R^2),
    title = paste("Enhanced Box Plot of R^2 for", trait),
    file_name = paste0("enhanced_boxplot_", trait, "_R2.png")
  )
  
  # Create Correlation plot
  create_enhanced_plot(
    data = trait_data,
    y_variable = "Correlation",
    y_label = "Correlation",
    title = paste("Enhanced Box Plot of Correlation for", trait),
    file_name = paste0("enhanced_boxplot_", trait, "_Correlation.png")
  )
}

print("Enhanced plots created and saved successfully!")
