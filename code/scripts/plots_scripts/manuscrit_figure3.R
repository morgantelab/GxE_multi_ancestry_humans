rm(list = ls())

# Load necessary libraries
library(ggplot2)
library(readr)
library(data.table)
library(dplyr)

# Define plot save directory
plot_save_dir <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/plots/"

# --- Load Eigen Data ---
eigenvec_data_pcrelate <- readRDS("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/pca_for_pcrelate.rds")

# Load phenotype/ethnicity data
load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/data2_20240930.RData")

# Merge with ethnicity data
eigenvec_pcrelate_df <- data.frame(ID = rownames(eigenvec_data_pcrelate$vectors),
                                   eigenvec_data_pcrelate$vectors,
                                   stringsAsFactors = FALSE)

# Ensure IDs are characters
eigenvec_pcrelate_df$ID <- as.character(eigenvec_pcrelate_df$ID)
dtt$ID <- as.character(dtt$ID)

# Merge datasets by ID
merged_data_pcrelate <- merge(eigenvec_pcrelate_df, dtt, by = "ID")

# Rename PCs clearly
colnames(merged_data_pcrelate)[2:ncol(eigenvec_pcrelate_df)] <- paste0("PC", 1:(ncol(eigenvec_pcrelate_df) - 1))

# Define ethnicity colors as specified
ethnicity_colors <- c(
  "White" = "limegreen",
          "Black" = "orange",
          "Asian" = "yellow",
          "Mixed" = "pink",
          "Chinese" = "red"
)

# Plot PC1 vs PC2 without title, formatted as manuscript
pc1_vs_pc2_plot <- ggplot(merged_data_pcrelate, aes(x = PC1, y = PC2, color = ethn1_consolidated)) +
  geom_point(alpha = 0.7, size = 1.5) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black"),
    axis.title = element_text(),
    legend.title = element_blank(),
    legend.position = "right"
  ) +
  scale_color_manual(values = ethnicity_colors) +
  labs(x = "Principal Component 1",
       y = "Principal Component 2")

# Save the plot
ggsave(filename = paste0(plot_save_dir, "Figure3_PC1_vs_PC2.png"),
       plot = pc1_vs_pc2_plot,
       width = 8, height = 6, dpi = 300)
