# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readr)
library(gridExtra)
library(grid)
# Load dataset
file_path <- "C:/Users/Aishwarya/Downloads/subset_dataset_figure1.csv"  
df <- read_csv(file_path)
# Convert ethnicity values to lowercase for consistency
df$ethn1_consolidated <- tolower(df$ethn1_consolidated)
# Define ethnicity color mapping
ethnicity_colors <- c(
  "white" = "limegreen",
          "black" = "orange",
          "asian" = "yellow",
          "mixed" = "pink",
          "chinese" = "red"
)
# Define common theme adjustments
common_theme <- theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10)  # Added margin to prevent cropping
  )
# Ensure y-axis scales are the same
y_limits <- range(c(df$DP0a, df$SP0a, df$PP0a), na.rm = TRUE)
# Define the order of ethnicity categories
ethnicity_order <- c("asian", "white", "black", "mixed", "chinese")
# Create figure (A) - Diastolic Blood Pressure (DP0a)
fig_a <- ggplot(df, aes(x = factor(ethn1_consolidated, levels = ethnicity_order), y = DP0a, fill = ethn1_consolidated)) +
  geom_boxplot(color = "black", alpha = 0.7) +
  scale_fill_manual(values = ethnicity_colors) +
  common_theme +
  labs(y = "Diastolic Blood Pressure", x = "") +  
  coord_cartesian(ylim = y_limits)
# Create figure (B) - Systolic Blood Pressure (SP0a)
fig_b <- ggplot(df, aes(x = factor(ethn1_consolidated, levels = ethnicity_order), y = SP0a, fill = ethn1_consolidated)) +
  geom_boxplot(color = "black", alpha = 0.7) +
  scale_fill_manual(values = ethnicity_colors) +
  common_theme +
  labs(y = "Systolic Blood Pressure", x = "") +  
  coord_cartesian(ylim = y_limits)
# Create figure (C) - Pulse Pressure (PP0a)
fig_c <- ggplot(df, aes(x = factor(ethn1_consolidated, levels = ethnicity_order), y = PP0a, fill = ethn1_consolidated)) +
  geom_boxplot(color = "black", alpha = 0.7) +
  scale_fill_manual(values = ethnicity_colors) +
  common_theme +
  labs(y = "Pulse Pressure", x = "") +  
  coord_cartesian(ylim = y_limits)
# Create figure (D) - Ethnicity Count Bar Plot
ethnicity_counts <- df %>%
  group_by(ethn1_consolidated) %>%
  summarise(count = n())
fig_d <- ggplot(ethnicity_counts, aes(x = factor(ethn1_consolidated, levels = ethnicity_order), y = count, fill = ethn1_consolidated)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = ethnicity_colors) +
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
ggsave("C:/Users/Aishwarya/Desktop/khushi didi/output/all_figures.png", final_plot, width = 8, height = 8, dpi = 300)
# Save all figures in one PDF with **STRICT SIZE**
pdf("C:/Users/Aishwarya/Desktop/khushi didi/output/all_figures.pdf", width = 8, height = 8)
grid.arrange(
  add_label(fig_a, "A"),
  add_label(fig_b, "B"),
  add_label(fig_c, "C"),
  add_label(fig_d, "D"),
  ncol = 2, nrow = 2
)
dev.off()