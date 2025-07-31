# Clean environment
rm(list=ls())

set.seed(1123)

# Load necessary libraries
library(ggplot2)
library(grid)
library(gridExtra)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(ggpubr)
library(cowplot)

# Load your CSV file
df <- read.csv("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/figure_4_dataset.csv", stringsAsFactors = FALSE)

# Map Component to friendly legend labels
component_labels <- c(
  "V_X1" = "Demographics",
  "V_X2" = "Structure",
  "V_G" = "Genetics",
  "V_E" = "Lifestyle",
  "V_GxE" = "Genetics-by-lifestyle"
)

# Keep only required components
df <- df %>% filter(Component %in% names(component_labels))

# Add readable Component label
df$Component_label <- component_labels[df$Component]

# Model mapping: Model number and component set for assignment
model_map <- c(
  "V_X1" = "Model 1: X1",
  "V_X1_X2" = "Model 2: X1, X2",
  "V_X1_X2_G" = "Model 3: X1, X2, G",
  "V_X1_X2_E" = "Model 4: X1, X2, E",
  "V_X1_X2_G_E" = "Model 5: X1, X2, G, E",
  "V_X1_X2_G_E_GxE" = "Model 6: X1, X2, G, E, GxE"
)

# Add new Model_full column for assignment
df$Model_full <- model_map[df$Model]

# For plotting, create a simple Model label for the x-axis
model_simple_map <- c(
  "Model 1: X1" = "Model 1",
  "Model 2: X1, X2" = "Model 2",
  "Model 3: X1, X2, G" = "Model 3",
  "Model 4: X1, X2, E" = "Model 4",
  "Model 5: X1, X2, G, E" = "Model 5",
  "Model 6: X1, X2, G, E, GxE" = "Model 6"
)
df$Model_simple <- model_simple_map[df$Model_full]

# Set model order for plotting
desired_model_full <- names(model_simple_map)
desired_model_simple <- unname(model_simple_map)

df$Model_full <- factor(df$Model_full, levels = desired_model_full)
df$Model_simple <- factor(df$Model_simple, levels = desired_model_simple)

# Ensure Component_label is ordered properly
df$Component_label <- factor(df$Component_label, levels = component_labels)

# Color palette for components (Genetics and Lifestyle colors swapped)
color_blind_palette <- c(
  "Demographics" = "#E69F00",
  "Structure" = "#56B4E9",
  "Genetics" = "purple",  # Swapped color
  "Lifestyle" = "#009E73", # Swapped color
  "Genetics-by-lifestyle" = "brown"
)

# Remove NA or zero values
df <- df %>% filter(!is.na(mean) & mean != 0)

# Flexible Y-limit function
get_ylim <- function(trait_name) {
  vals <- df %>% filter(trait == trait_name) %>% pull(ci_upper)
  if(length(vals) == 0) return(c(0, 1))
  max_val <- max(vals, na.rm = TRUE)
  c(0, max_val * 1.05)
}

# Trait plotting function with optional custom ylim
plot_trait <- function(trait_name, legend_text_size=35, custom_ylim=NULL) {
  ylim <- if(is.null(custom_ylim)) get_ylim(trait_name) else custom_ylim
  ggplot(df %>% filter(trait == trait_name),
         aes(x = Model_simple, y = mean, color = Component_label)) +
    geom_point(position = position_dodge2(width = 1, preserve = "single"), size = 9) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 1,
                  position = position_dodge2(width = 1, preserve = "single")) +
    labs(title = "", x = "", y = "Variance") +
    scale_y_continuous(limits = ylim, expand = c(0, 0)) +
    scale_x_discrete(labels = desired_model_simple, expand = expansion(add = 0.6)) +
    scale_color_manual(values = color_blind_palette,
                       labels = names(color_blind_palette)) +
    theme_bw(base_size = 15) +  # White background
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = legend_text_size),
      legend.key.size = unit(2, "cm"),
      legend.background = element_rect(fill = "white", color = NA),
      legend.box.background = element_rect(fill = "white", color = NA),
      axis.text.x = element_text(size = 35, color = "black", angle = 0, hjust = 0.5),
      axis.text.y = element_text(size = 35, color = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 35, margin = margin(t = 0, r = 20, b = 0, l = 0)),
      plot.margin = unit(c(1.5, 0.8, 0.8, 0.8), "cm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.ticks.length = unit(0.3, "cm"),
      panel.border = element_rect(colour = "black", fill=NA, size=2),
      plot.background = element_rect(fill = "white", color = NA)
    )
}

# Generate individual trait plots
plot_dp <- plot_trait("DP", legend_text_size=35)
plot_sp <- plot_trait("SP", legend_text_size=30)
plot_pp <- plot_trait("PP", legend_text_size=30, custom_ylim=c(-0.5, get_ylim("PP")[2]))

# Arrange all plots vertically with common legend
plots_arranged <- ggarrange(plot_dp, plot_sp, plot_pp, ncol=1,
                            common.legend = TRUE, legend = "bottom")

# Add A, B, C labels to the left
final <- ggdraw(plots_arranged) +
  draw_plot_label(
    label = c("A", "B", "C"),
    x = rep(0.01, 3),
    y = c(0.98, 0.65, 0.32),
    hjust = 0, vjust = 1,
    fontface = "bold", size = 50, color = "black"
  )

# Save the final PNG with white background
ggsave("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/plots/Fig4_var_exp_ci.pdf", final, height = 32, width = 24, dpi = 300, device = "pdf")
