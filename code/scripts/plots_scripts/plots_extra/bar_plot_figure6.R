rm(list=ls()); gc()
set.seed(1123)

# Load required libraries
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(purrr)
library(cowplot)

# Define paths
summary_file <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/all_model_summary_pval_ethn_corrected.csv"

# Load summary data
summary_data <- read_csv(summary_file)

# Filter out rows with GEselect in the Model name
summary_data <- summary_data %>%
  filter(!grepl("GEselect", Model))

# Define traits and ancestries for ordering
traits <- c("SP", "DP", "PP")
ancestries <- c("asian", "mixed", "black", "chinese", "white")
models <- c("X1", "X1_X2", "X1_X2_G", "X1_X2_E", "X1_X2_G_E", "X1_X2_G_E_GE")

# Format data
summary_data_ancestries <- summary_data %>%
  filter(Trait %in% traits, Model %in% models) %>%
  mutate(
    Trait = factor(Trait, levels = traits),
    Ancestry = factor(Ancestry, levels = ancestries),
    Model = factor(Model, levels = models)
  )

# Define ethnicity colors
ethnicity_colors <- c(
  "white" = "limegreen",
          "black" = "orange",
          "asian" = "yellow",
          "mixed" = "pink",
          "chinese" = "red"
)

# Create labeled trait titles
traits_titles <- c("DP" = "A", "SP" = "B", "PP" = "C")

# Find global y-axis limit
y_max <- max(summary_data_ancestries$R_squared, na.rm = TRUE)

# Create bar plots with manuscript styling and consistent y-axis
plot_list <- lapply(names(traits_titles), function(trait) {
  ggplot(summary_data_ancestries %>% filter(Trait == trait), aes(x = Model, y = R_squared, fill = Ancestry)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black", width = 0.9) +
    scale_fill_manual(values = ethnicity_colors) +
    theme_bw(base_size = 16) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(color = "black", size = 14),
      axis.title.y = element_text(face = "bold", size = 16),
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 30, hjust = 1),
      legend.position = "none"
    ) +
    labs(
      title = traits_titles[[trait]],
      y = expression(R^2)
    ) +
    ylim(0, y_max)
})

# Create one legend from one plot
legend <- get_legend(
  plot_list[[1]] + theme(legend.position = "bottom", legend.title = element_blank())
)

# Combine plots into one column
combined_plot <- plot_grid(plotlist = plot_list, ncol = 1, align = "v")
final_plot <- plot_grid(combined_plot, legend, ncol = 1, rel_heights = c(1, 0.1))

# Save combined figure with increased size
ggsave("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/plots/barplot_preds_models_ancestries_combined_column.png",
       plot = final_plot,
       width = 12, height = 18, dpi = 300)
