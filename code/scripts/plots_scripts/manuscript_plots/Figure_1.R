rm(list=ls()); gc()
set.seed(1123)

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readr)
library(gridExtra)
library(grid)

# Load dataset
file_path <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/subset_dataset_figure1.csv"
df <- read_csv(file_path)

# Capitalize only the first letter of each ethnicity value
df$ethn1_consolidated <- paste0(toupper(substring(df$ethn1_consolidated, 1, 1)),
                                tolower(substring(df$ethn1_consolidated, 2)))

# Define ethnicity color mapping
ethnicity_colors <- c(
  "White" = "limegreen",
  "Black" = "orange",
  "Asian" = "yellow",
  "Mixed" = "pink",
  "Chinese" = "red"
)

# Define common theme adjustments (with plot box border)
common_theme <- theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10)  # Added margin to prevent cropping
  )

# Define the order of ethnicity categories
ethnicity_order <- c("Asian", "White", "Black", "Mixed", "Chinese")

# Create figure (A) - Diastolic Blood Pressure (DP0a)
fig_a <- ggplot(df, aes(x = factor(ethn1_consolidated, levels = ethnicity_order), y = DP0a, fill = ethn1_consolidated)) +
  geom_boxplot(color = "black", alpha = 0.7) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3.5,
               fill = "deepskyblue3", color = "black") +
  scale_fill_manual(values = ethnicity_colors) +
  scale_x_discrete(labels = ethnicity_order) +
  common_theme +
  labs(y = "mmHg", x = "")

# Create figure (B) - Systolic Blood Pressure (SP0a)
fig_b <- ggplot(df, aes(x = factor(ethn1_consolidated, levels = ethnicity_order), y = SP0a, fill = ethn1_consolidated)) +
  geom_boxplot(color = "black", alpha = 0.7) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3.5,
               fill = "deepskyblue3", color = "black") +
  scale_fill_manual(values = ethnicity_colors) +
  scale_x_discrete(labels = ethnicity_order) +
  common_theme +
  labs(y = "mmHg", x = "")

# Create figure (C) - Pulse Pressure (PP0a) with improved y-axis
fig_c <- ggplot(df, aes(x = factor(ethn1_consolidated, levels = ethnicity_order), y = PP0a, fill = ethn1_consolidated)) +
  geom_boxplot(color = "black", alpha = 0.7) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3.5,
               fill = "deepskyblue3", color = "black") +
  scale_fill_manual(values = ethnicity_colors) +
  scale_x_discrete(labels = ethnicity_order) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +  # improved spacing
  common_theme +
  labs(y = "mmHg", x = "")

# Create figure (D) - Ethnicity Count Bar Plot (with flexible y-axis)
ethnicity_counts <- df %>%
  group_by(ethn1_consolidated) %>%
  summarise(count = n())

fig_d <- ggplot(ethnicity_counts, aes(x = factor(ethn1_consolidated, levels = ethnicity_order), y = count, fill = ethn1_consolidated)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = ethnicity_colors) +
  scale_x_discrete(labels = ethnicity_order) +
  common_theme +
  labs(y = "Count", x = "")

# Function to add labels (A, B, C, D) **INSIDE** the plot
add_label <- function(plot, label) {
  arrangeGrob(
    textGrob(label, x = 0.1, y = 0.50, just = "left", gp = gpar(fontsize = 14, fontface = "bold")),
    plot,
    heights = c(0.05, 1),  # Reduced label height
    ncol = 1
  )
}

# Arrange figures with properly positioned labels **INSIDE IMAGE BOUNDS**
final_plot <- grid.arrange(
  add_label(fig_a, "A"),
  add_label(fig_b, "B"),
  add_label(fig_c, "C"),
  add_label(fig_d, "D"),
  ncol = 2, nrow = 2
)

# Save all figures as PNG with **STRICT SIZE**
#ggsave("C:/Users/Aishwarya/Desktop/khushi didi.png", final_plot, width = 8, height = 8, dpi = 300)

# Save all figures in one PDF with **STRICT SIZE**
pdf ("C:/Users/Aishwarya/Desktop/all.pdf", width = 8, height = 8)
grid.arrange(
  add_label(fig_a, "A"),
  add_label(fig_b, "B"),
  add_label(fig_c, "C"),
  add_label(fig_d, "D"),
  ncol = 2, nrow = 2
)
dev.off()
