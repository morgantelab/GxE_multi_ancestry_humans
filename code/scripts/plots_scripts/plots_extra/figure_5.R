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
summary_file <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/all_model_summary_folds_pval.csv"

# Load summary data
summary_data <- read_csv(summary_file)

# Filter out GEselect rows if present
summary_data <- summary_data %>% filter(!grepl("GEselect", Model))

# Define parameters
traits <- c("SP", "DP", "PP")
models <- c("X1", "X1_X2", "X1_X2_G", "X1_X2_E", "X1_X2_G_E", "X1_X2_G_E_GE")
model_labels <- c("X1" = "Model 1", "X1_X2" = "Model 2", "X1_X2_G" = "Model 3", 
                  "X1_X2_E" = "Model 4", "X1_X2_G_E" = "Model 5", "X1_X2_G_E_GE" = "Model 6")
folds <- 1:5

# Format data
summary_data_folds <- summary_data %>%
  filter(Trait %in% traits, Model %in% models) %>%
  mutate(
    Trait = factor(Trait, levels = traits),
    Fold = factor(Fold, levels = folds),
    Model = factor(Model, levels = models)
  )

# Define trait colors
trait_colors <- c("DP" = "limegreen", "SP" = "blue", "PP" = "orange")
traits_titles <- c("DP" = "A", "SP" = "B", "PP" = "C")

# Determine global y-axis limit
y_max <- max(summary_data_folds$R_squared, na.rm = TRUE)

# Create updated boxplots with custom model labels
plot_list <- lapply(names(traits_titles), function(trait) {
  ggplot(summary_data_folds %>% filter(Trait == trait),
         aes(x = Model, y = R_squared)) +
    geom_boxplot(color = "black", alpha = 0.6, fill = trait_colors[trait], width = 0.6) +
    scale_x_discrete(labels = model_labels) +
    theme_bw(base_size = 28) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(color = "black", size = 28),
      axis.title.x = element_blank(),
      axis.title.y = element_text(face = "bold", size = 30),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, margin = margin(t = 12)),
      legend.position = "none"
    ) +
    labs(
      title = traits_titles[[trait]],
      y = expression(R^2)
    ) +
    ylim(0, y_max)
})

# Combine plots into one column
final_plot <- plot_grid(plotlist = plot_list, ncol = 1, align = "v")

# Save final plot
ggsave("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/plots/boxplot_preds_models_folds_combined_column.png",
       plot = final_plot,
       width = 16, height = 24, dpi = 300)
