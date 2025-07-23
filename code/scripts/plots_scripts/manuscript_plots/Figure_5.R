rm(list=ls()); gc()
set.seed(1123)

# Load required libraries
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(purrr)
library(cowplot)
library(scales)

# Load the summary file
summary_file <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/all_model_summary_ethn_corrected.csv"
summary_data <- read_csv(summary_file)

# Define levels
traits <- c("SP", "DP", "PP")
ancestries <- c("asian", "mixed", "black", "chinese", "white")
model_levels <- c("X1", "X1_X2", "X1_X2_G", "X1_X2_E", "X1_X2_G_E", "X1_X2_G_E_GE")
model_labels <- c("Model 1", "Model 2", "Model 3", "Model 4", "Model 5", "Model 6")

# Filter and format data
summary_data_ancestries <- summary_data %>%
  filter(Trait %in% traits, Model %in% model_levels) %>%
  mutate(
    Trait = factor(Trait, levels = traits),
    Ancestry = factor(str_to_title(Ancestry), levels = str_to_title(ancestries)),
    Model = factor(Model, levels = model_levels),
    models = factor(models, levels = model_labels)
  )

# Define custom model colors to match the image
model_colors <- c(
  "X1" = "#1b9e77",        # Teal/Dark Green (Model 1)
  "X1_X2" = "#d95f02",       # Orange (Model 2)
  "X1_X2_G" = "#7570b3",    # Purple (Model 3)
  "X1_X2_E" = "#e7298a",    # Pink/Magenta (Model 4)
  "X1_X2_G_E" = "#66a61e",   # Lighter Green (Model 5) - adjusted to be distinct from M1, based on the image
  "X1_X2_G_E_GE" = "#e6ab02" # Yellow/Gold (Model 6)
)

# Trait plot labels
traits_titles <- c("DP" = "A", "SP" = "B", "PP" = "C")

# Create individual bar plots with flexible y-axis
plot_list <- lapply(names(traits_titles), function(trait) {
  ggplot(summary_data_ancestries %>% filter(Trait == trait),
         aes(x = Ancestry, y = R_squared, fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8),
             color = "black", width = 0.7) +
    scale_fill_manual(
      values = model_colors,
      breaks = model_levels,
      labels = model_labels
    ) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
    theme_bw(base_size = 20) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(color = "black", size = 40), # Increased size of x-axis labels
      axis.text.y = element_text(color = "black", size = 40), # Increased size of y-axis labels
      axis.title.y = element_text(face = "bold", size = 40), # Increased size of y-axis title to match labels
      axis.title.x = element_blank(),
      legend.position = "none",
      plot.title.position = "plot",
      plot.title = element_text(hjust = 0, face = "bold", size = 60, color = "black"), # Adjusted hjust to 0 for more left alignment
      panel.border = element_rect(colour = "black", fill=NA, size=3) # Increased size of box outline
    ) +
    labs(
      title = traits_titles[[trait]],
      y = expression(italic(R^2))
    )
})

# Stack the plots vertically
combined_plots <- plot_grid(plotlist = plot_list, ncol = 1, align = "v")

# Create a dummy plot just for the legend
legend_plot <- ggplot(summary_data_ancestries, aes(x = Ancestry, y = R_squared, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           color = "black", width = 0.7) +
  scale_fill_manual(
    values = model_colors,
    breaks = model_levels,
    labels = model_labels
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 40), # Increased size of legend text
    legend.key.size = unit(1.5, "cm")
  )

# Extract the legend
legend <- get_legend(legend_plot)

# Combine the three plots vertically
all_plots <- plot_grid(plot_list[[1]], plot_list[[2]], plot_list[[3]],
                       ncol = 1, align = "v", label_size = 20)

# Combine plots with legend side by side
final_plot <- plot_grid(all_plots, legend, rel_widths = c(0.9, 0.15))

# Save the final PDF
ggsave("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/plots/Figure_5.pdf", height = 30, width = 24, dpi = 300, device = "pdf")