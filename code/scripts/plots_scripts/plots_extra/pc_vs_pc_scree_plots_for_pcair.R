# Load the required data
load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/data2_20240930.RData")
new_dataset <- dtt
pcair_results <- readRDS("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/pcair_for_pcrelate_results.rds")

# Extract the first three principal components from pcair_results
pca_scores <- data.frame(SampleID = pcair_results$sample.id, 
                         PC1 = pcair_results$vectors[, 1], 
                         PC2 = pcair_results$vectors[, 2], 
                         PC3 = pcair_results$vectors[, 3])

# Merge PCA scores with ethnicity data from new_dataset
combined_data <- merge(pca_scores, new_dataset, by.x = "SampleID", by.y = "ID")

# Define the font sizes
title_size <- 18
axis_title_size <- 16
axis_text_size <- 14

# Plot PC1 vs PC2
plot_pc1_vs_pc2 <- ggplot(combined_data, aes(x = PC1, y = PC2, color = ethn1_consolidated)) +
  geom_point(alpha = 0.6) +
  labs(title = "PC1 vs PC2 Colored by Ethnicity", x = "PC1", y = "PC2") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = title_size),
    axis.title = element_text(size = axis_title_size),
    axis.text = element_text(size = axis_text_size)
  )

# Save the plot
ggsave(filename = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/plots/PC1_vs_PC2_ethn_pcair.png",
       plot = plot_pc1_vs_pc2, width = 10, height = 8, units = "in")

# Plot PC2 vs PC3
plot_pc2_vs_pc3 <- ggplot(combined_data, aes(x = PC2, y = PC3, color = ethn1_consolidated)) +
  geom_point(alpha = 0.6) +
  labs(title = "PC2 vs PC3 Colored by Ethnicity", x = "PC2", y = "PC3") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = title_size),
    axis.title = element_text(size = axis_title_size),
    axis.text = element_text(size = axis_text_size)
  )

# Save the plot
ggsave(filename = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/plots/PC2_vs_PC3_ethn_pcair.png",
       plot = plot_pc2_vs_pc3, width = 10, height = 8, units = "in")

# Plot PC1 vs PC3
plot_pc1_vs_pc3 <- ggplot(combined_data, aes(x = PC1, y = PC3, color = ethn1_consolidated)) +
  geom_point(alpha = 0.6) +
  labs(title = "PC1 vs PC3 Colored by Ethnicity", x = "PC1", y = "PC3") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = title_size),
    axis.title = element_text(size = axis_title_size),
    axis.text = element_text(size = axis_text_size)
  )

# Save the plot
ggsave(filename = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/plots/PC1_vs_PC3_ethn_pcair.png",
       plot = plot_pc1_vs_pc3, width = 10, height = 8, units = "in")

# Create a data frame for the scree plot
var_explained <- pcair_results$values
proportion_explained <- var_explained / sum(var_explained) * 100
cumulative_explained <- cumsum(proportion_explained)

pc_data <- data.frame(PC = 1:length(proportion_explained), 
                      Variance = proportion_explained, 
                      CumulativeVariance = cumulative_explained)

# Standard Scree Plot
standard_scree_plot <- ggplot(pc_data, aes(x = PC, y = Variance)) +
  geom_point() +  # Add points
  geom_line() +  # Connect the points with lines
  labs(title = "Scree Plot", x = "Principal Component", y = "Percentage of Variance Explained") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = title_size),
    axis.title = element_text(size = axis_title_size),
    axis.text = element_text(size = axis_text_size)
  )

# Save the Standard Scree Plot
ggsave(filename = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/plots/Standard_Scree_Plot_pcair.png",
       plot = standard_scree_plot, width = 10, height = 8, units = "in")

# Cumulative Scree Plot
cumulative_scree_plot <- ggplot(pc_data, aes(x = PC, y = CumulativeVariance)) +
  geom_point(color = "red") +  # Add points
  geom_line(color = "red") +  # Connect the points with lines
  labs(title = "Cumulative Scree Plot", x = "Principal Component", y = "Cumulative Percentage of Variance Explained") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = title_size),
    axis.title = element_text(size = axis_title_size),
    axis.text = element_text(size = axis_text_size)
  )

# Save the Cumulative Scree Plot
ggsave(filename = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/plots/Cumulative_Scree_Plot_pcair.png",
       plot = cumulative_scree_plot, width = 10, height = 8, units = "in")
