rm(list = ls())

# Load necessary libraries
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(cowplot)

# # Load the variance components data
# variance_data <- read_csv("/mnt/data/Variance_Components_Summary.csv")
Variance_Components_Summary <- read_csv("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/Variance_Components_Summary.csv")


# Add new X1_X2_E results (Model 4) manually for DP, SP, PP
new_data <- tibble::tribble(
  ~MODEL, ~filename, ~V_E, ~V_G, ~V_GxE, ~V_X1, ~V_X2,
  "X1_X2_E", "DP", 10.4319042, NA, NA, 4.1769837, 1.4419572,
  "X1_X2_E", "SP", 5.0994285,   NA, NA, 13.8046433, 0.8159450,
  "X1_X2_E", "PP", 1.1946591,   NA, NA, 19.6632455, 0.3320431
)

# Merge new data into the original dataset
variance_data <- bind_rows(Variance_Components_Summary, new_data)

# Transform data into long format for ggplot
variance_long <- variance_data %>%
  pivot_longer(cols = c(V_E, V_G, V_GxE, V_X1, V_X2),
               names_to = "Component",
               values_to = "Variance") %>%
  filter(!is.na(Variance))

# Rename models for clarity
model_labels <- c(
  "just_X1" = "Model 1",
  "just_X1_X2" = "Model 2",
  "X1_X2_G" = "Model 3",
  "X1_X2_E" = "Model 4",
  "X1_X2_G_E" = "Model 5",
  "X1_X2_G_E_GE" = "Model 6"
)
variance_long$MODEL <- factor(variance_long$MODEL, levels = names(model_labels), labels = model_labels)

# Set desired order and colorblind-friendly colors for components
component_order <- c("V_X1", "V_X2", "V_E", "V_G", "V_GxE")
component_labels <- c("X1", "X2", "E", "G", "GxE")
component_colors <- c(
  "X1" = "#D55E00",
  "X2" = "#E69F00",
  "E"  = "#009E73",
  "G"  = "#0072B2",
  "GxE"= "#CC79A7"
)

variance_long$Component <- factor(variance_long$Component, levels = component_order, labels = component_labels)

# Trait titles
traits_titles <- list("DP" = "A", "SP" = "B", "PP" = "C")
plots_list <- list()

for (trait in names(traits_titles)) {
  trait_data <- variance_long %>% filter(filename == trait)
  
  base_plot <- ggplot(trait_data, aes(x = MODEL, y = Variance, color = Component)) +
    geom_point(size = 4, position = position_dodge(width = 0.5)) +
    scale_color_manual(values = component_colors, drop = FALSE) +
    theme_bw(base_size = 16) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      legend.position = "none",
      aspect.ratio = 1
    ) +
    labs(x = "Model",
         y = "Variance",
         title = traits_titles[[trait]])
  
  plots_list[[trait]] <- base_plot
}

# Create clean dummy plot for extracting legend
legend_plot <- ggplot(variance_long, aes(x = MODEL, y = Variance, color = Component)) +
  geom_point(size = 4) +
  scale_color_manual(values = component_colors, drop = FALSE) +
  theme_void(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.key.width = unit(2, "lines"),
    legend.key.height = unit(1, "lines")
  )

legend_grob <- get_legend(legend_plot)

# Combine row and legend separately
combined_row <- plot_grid(plots_list$DP, plots_list$SP, plots_list$PP, nrow = 1, align = "hv", labels = NULL)
final_plot <- plot_grid(combined_row, legend_grob, ncol = 1, rel_heights = c(1, 0.1))

# Save final plot
ggsave(filename = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/plots/Variance_Components_Grid_Row_legend_withE.png",
       plot = final_plot,
       width = 18, height = 8, dpi = 300)
