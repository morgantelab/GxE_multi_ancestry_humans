rm(list = ls())

# Load necessary libraries
library(ggplot2)
library(readr)
library(data.table)
library(dplyr)
library(cowplot)

# Define plot save directory
plot_save_dir <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/plots/"

# Load phenotype/ethnicity data
dtt <- readRDS("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/data3_20250514.rds")
dtt$ID <- as.character(dtt$ID)

# Define ethnicity colors
ethnicity_colors <- c(
  "White" = "limegreen",
  "Black" = "orange",
  "Asian" = "yellow",
  "Mixed" = "pink",
  "Chinese" = "red"
)

# Shared base theme with bigger text
base_theme <- theme_bw(base_size = 50) +  # Bigger overall text
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black", size = 40),
    axis.title = element_text(size = 42),
    legend.text = element_text(size = 40),
    plot.title = element_text(size = 40)
  )

# --- Panel A: PLINK PCA ---
eigenvec_data_plink <- readRDS("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/pca_for_plink.rds")
eigenvec_plink_df <- data.frame(ID = rownames(eigenvec_data_plink$vectors),
                                eigenvec_data_plink$vectors,
                                stringsAsFactors = FALSE)
eigenvec_plink_df$ID <- as.character(eigenvec_plink_df$ID)
merged_data_plink <- merge(eigenvec_plink_df, dtt, by = "ID")
colnames(merged_data_plink)[2:ncol(eigenvec_plink_df)] <- paste0("PC", 1:(ncol(eigenvec_plink_df) - 1))

plot_plink <- ggplot(merged_data_plink, aes(x = PC1, y = PC2, color = ethn1_consolidated)) +
  geom_point(alpha = 1, size = 4) +
  base_theme +
  theme(legend.position = "none") +
  scale_color_manual(values = ethnicity_colors) +
  labs(x = "Principal Component 1", y = "Principal Component 2")

# --- Panel B: PC-Relate PCA ---
eigenvec_data_pcrelate <- readRDS("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/pca_for_pcrelate.rds")
eigenvec_pcrelate_df <- data.frame(ID = rownames(eigenvec_data_pcrelate$vectors),
                                   eigenvec_data_pcrelate$vectors,
                                   stringsAsFactors = FALSE)
eigenvec_pcrelate_df$ID <- as.character(eigenvec_pcrelate_df$ID)
merged_data_pcrelate <- merge(eigenvec_pcrelate_df, dtt, by = "ID")
colnames(merged_data_pcrelate)[2:ncol(eigenvec_pcrelate_df)] <- paste0("PC", 1:(ncol(eigenvec_pcrelate_df) - 1))

plot_pcrelate_no_legend <- ggplot(merged_data_pcrelate, aes(x = PC1, y = PC2, color = ethn1_consolidated)) +
  geom_point(alpha = 1, size = 4) +
  base_theme +
  theme(legend.position = "none") +
  scale_color_manual(values = ethnicity_colors) +
  labs(x = "Principal Component 1", y = "Principal Component 2")

# Extract legend separately
# Extract legend separately with better spacing
plot_for_legend <- ggplot(merged_data_pcrelate, aes(x = PC1, y = PC2, color = ethn1_consolidated)) +
  geom_point(size = 6) +
  theme_minimal(base_size = 50) +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 40),
    legend.key.height = unit(4, "lines")  # ⬅️ more vertical spacing
  ) +
  scale_color_manual(
    values = ethnicity_colors,
    guide = guide_legend(
      override.aes = list(size = 4),  # ⬅️ larger points in legend
      keyheight = unit(4, "lines")  # ⬅️ match vertical space
    )
  )


legend_shared <- get_legend(plot_for_legend)

# Stack plots vertically with labels A and B
stacked_plots <- plot_grid(
  plot_plink, plot_pcrelate_no_legend,
  labels = c("A", "B"),
  label_size = 55,
  ncol = 1,
  align = "v"
)

# Combine stacked plots + legend to the right
final_plot <- plot_grid(
  stacked_plots, legend_shared,
  ncol = 2,
  rel_widths = c(1, 0.3)
)

# Save final figure
ggsave(filename = file.path(plot_save_dir, "Figure2_PCA_PLINK_vs_PCRELATE_vertical_sharedlegend_biggertext.pdf"),
       plot = final_plot,
       width = 24, height = 32)
