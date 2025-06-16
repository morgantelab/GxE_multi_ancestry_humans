# Clean environment
rm(list=ls())
# Load necessary libraries
library(ggplot2)
library(grid)
library(gridExtra)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(ggpubr)
library(cowplot)  # For placing labels outside the plot

# Load your CSV file
df <- read.csv("/data2/morgante_lab/ukbiobank_projects/kgoda_practice/gwas_analysis/cohort_batch_effects/plots_manuscript/sep_run_2024/VCEval_20240823.xlsx", stringsAsFactors = FALSE)

# Map Component to friendly legend labels
component_labels <- c(
  "V_X1" = "Demographics",
  "V_X2" = "Structure",
  "V_G" = "Genetics",
  "V_E" = "Lifestyle",
  "V_GxE" = "Genetics-by-lifestyle"
)

# Only keep rows with components of interest
df <- df %>% filter(Component %in% names(component_labels))

# Add a pretty label column for plotting
df$Component_label <- component_labels[df$Component]

# Recode Model to Model 1, Model 2, ... in order of appearance
model_levels <- unique(df$Model)
model_labels <- paste0("Model ", seq_along(model_levels))
model_map <- setNames(model_labels, model_levels)
df$Model <- model_map[df$Model]

# Set the desired order for the x-axis (Model 1 to Model 5)
desired_model_labels <- paste0("Model ", 1:6)
df$Model <- factor(df$Model, levels = desired_model_labels)

# Ensure Component_label is a factor and ordered as in the legend
df$Component_label <- factor(df$Component_label, levels = component_labels)

# X-axis labels as "Model 1", "Model 2", ...
Lables <- desired_model_labels

# Color-blind-friendly palette (Okabe-Ito) for the components
color_blind_palette <- c(
  "Demographics" = "#E69F00",        # orange
  "Structure" = "#56B4E9",          # sky blue
  "Genetics" = "#009E73",           # green
  "Lifestyle" = "#F39B7F",          # reddish orange
  "Genetics-by-lifestyle" = "#0072B2" # blue
)

# Filter out rows with 0 or NA mean values
df <- df %>% filter(!is.na(Mean) & Mean != 0)

# Function to get flexible y-axis limits for each trait
get_ylim <- function(trait_name) {
  vals <- df %>% filter(Trait == trait_name) %>% pull(CI_Upper)
  if(length(vals) == 0) return(c(0, 1))
  max_val <- max(vals, na.rm = TRUE)
  c(0, max_val * 1.05)  # Slight buffer to prevent cutoff
}

# Create individual plots for each trait (no panel label here)
plot_trait <- function(trait_name, legend_text_size=35) {
  ylim <- get_ylim(trait_name)
  ggplot(df %>% filter(Trait == trait_name), 
         aes(x = Model, y = Mean, color = Component_label)) +
    geom_point(position = position_dodge2(width = 1, preserve = "single"), size = 9) +
    geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 1, 
                  position = position_dodge2(width = 1, preserve = "single")) +
    labs(title = "", x = "", y = "Variance") +
    scale_y_continuous(limits = ylim, expand = c(0, 0)) +
    scale_x_discrete(labels = Lables, expand = expansion(add = 0.6)) +
    scale_color_manual(values = color_blind_palette, 
                       labels = names(color_blind_palette)) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = legend_text_size),
      legend.key.size = unit(2, "cm"),
      axis.text.x = element_text(size = 35, color = "black"),
      axis.text.y = element_text(size = 35, color = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 35, margin = margin(t = 0, r = 20, b = 0, l = 0)),
      plot.margin = unit(c(1.5, 0.8, 0.8, 0.8), "cm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.ticks.length = unit(0.3, "cm"),
      panel.border = element_rect(colour = "black", fill=NA, size=2)
    )
}

plot_dp <- plot_trait("DP", legend_text_size=35)
plot_sp <- plot_trait("SP", legend_text_size=30)
plot_pp <- plot_trait("PP", legend_text_size=30)

# Arrange plots with common legend (no panel labels here)
plots_arranged <- ggarrange(plot_dp, plot_sp, plot_pp, ncol=1, 
                            common.legend = TRUE, legend = "bottom")

# Use cowplot to add A, B, C outside the left border
final <- ggdraw(plots_arranged) +
  draw_plot_label(
    label = c("A", "B", "C"),
    x = rep(0.01, 3), # far left
    y = c(0.98, 0.65, 0.32), # adjust for 3 panels
    hjust = 0, vjust = 1,
    fontface = "bold", size = 50, color = "black"
  )

ggsave("Fig2_var_exp_BP_noah_sugg.pdf", final, height = 32, width = 24, dpi = 300, device = "pdf")