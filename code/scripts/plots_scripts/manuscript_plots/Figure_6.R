# Clear environment
rm(list = ls()); gc()
set.seed(1123)

# Load required libraries
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(purrr)
library(cowplot)

# Load updated data
summary_file <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/all_model_summary_folds.csv"
summary_data <- read_csv(summary_file)

# Filter out GEselect rows
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


# Colorblind-friendly color palette (Set2)
trait_colors <- c("DP" = "#66C2A5", "SP" = "#8DA0CB", "PP" = "#FC8D62")
traits_titles <- c("DP" = "A", "SP" = "B", "PP" = "C")

# Generate boxplots with flexible y-axis per trait
plot_list <- lapply(names(traits_titles), function(trait) {
  data_trait <- summary_data_folds %>% filter(Trait == trait)
  
  # Calculate mean R² values for diamonds
  mean_data <- data_trait %>%
    group_by(Model) %>%
    summarize(mean_r2 = mean(R_squared, na.rm = TRUE), .groups = 'drop')
  
  p <- ggplot(data_trait, aes(x = Model, y = R_squared)) +
    geom_boxplot(color = "black", alpha = 0.6, fill = trait_colors[trait], width = 0.6) +
    geom_point(data = mean_data, aes(x = Model, y = mean_r2),
               shape = 23, size = 4.5, fill = "#FFD700", color = "black", stroke = 1.2) +  # Filled golden diamonds
    scale_x_discrete(labels = model_labels) +
    theme_bw(base_size = 28) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      axis.line = element_blank(),
      axis.text = element_text(color = "black", size = 28),
      axis.title.x = element_blank(),
      axis.title.y = element_text(face = "bold", size = 30),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, margin = margin(t = 12)),
      legend.position = "none",
      plot.title = element_text(hjust = -0.1, face = "bold")
    ) +
    labs(
      title = traits_titles[[trait]],
      y = expression(italic(R^2))
    )
  
  # Custom y-axis for PP and SP
  if (trait == "PP") {
    p <- p + scale_y_continuous(
      limits = c(0.19, 0.235),
      breaks = seq(0.19, 0.235, by = 0.01),
      labels = scales::number_format(accuracy = 0.01)
    )
  } else if (trait == "SP") {
    # Auto-determine useful breaks with gap of 0.03
    y_min <- floor(min(data_trait$R_squared, na.rm = TRUE) * 50) / 50  # closest lower multiple of 0.02
    y_max <- ceiling(max(data_trait$R_squared, na.rm = TRUE) * 50) / 50  # closest upper multiple of 0.02
    p <- p + scale_y_continuous(
      breaks = seq(y_min, y_max, by = 0.03),
      labels = scales::number_format(accuracy = 0.01)
    )
  } else {
    p <- p + scale_y_continuous(labels = scales::number_format(accuracy = 0.01))
  }
  
  return(p)
})

# Combine plots into a column
final_plot <- plot_grid(plotlist = plot_list, ncol = 1, align = "v")

# Save final plot
ggsave("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/plots/Figure_6_boxplot.pdf",
       plot = final_plot,
       width = 16, height = 24, dpi = 300)

