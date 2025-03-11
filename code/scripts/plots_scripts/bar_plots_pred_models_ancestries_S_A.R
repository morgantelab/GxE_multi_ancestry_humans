# Clear workspace and set seed for reproducibility
rm(list=ls()); gc()
set.seed(1123)

# Load required libraries
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(purrr)

# Define paths
summary_file <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/Updated_predictions_summary_S_A.csv"

# Load summary data
summary_data <- read_csv(summary_file)

# Define parameters used to generate filenames
traits <- c("SP", "DP", "PP")
models <- c("S_A_X2", "S_A_X2_G", "S_A_X2_E", "S_A_X2_G_E", "S_A_X2_G_E_GE")
ancestries <- c("asian", "mixed", "black", "chinese", "white")  # Replace folds with ancestries

# Generate the correct ancestry-based dataset with explicit Trait, Ancestry, and Model
ancestry_data <- expand.grid(Trait = traits, Ancestry = ancestries, Model = models, stringsAsFactors = FALSE) %>%
  mutate(File = sprintf("PREDs_%s_ethn_%s_%s.csv", Trait, Ancestry, Model))

# Merge with summary data to fill in R² values
summary_data_ancestries <- ancestry_data %>%
  left_join(summary_data, by = "File") %>%  # Merging to fill R² values
  select(Trait, Ancestry, Model, R_Squared)  # Keep only necessary columns

# Convert Ancestry to factor for proper ordering
summary_data_ancestries <- summary_data_ancestries %>%
  mutate(Ancestry = factor(Ancestry, levels = ancestries))

# Verify extraction correctness
print(table(summary_data_ancestries$Trait))  # Ensure Trait values are correct
print(table(summary_data_ancestries$Model))  # Ensure Model values are correct
print(table(summary_data_ancestries$Ancestry))   # Ensure Ancestry values are correct

# Define colors for each trait
trait_colors <- c("DP" = "limegreen", "SP" = "blue", "PP" = "orange")

# Create bar plots with styling
plot_list <- lapply(unique(summary_data_ancestries$Trait), function(trait) {
  ggplot(summary_data_ancestries %>% filter(Trait == trait), aes(x = Model, y = R_Squared, fill = Ancestry)) +
    geom_bar(stat = "identity", position = "dodge", color = "black") +  # Grouped bar plot with black outline
    scale_fill_manual(values = c("asian" = "lightblue", "mixed" = "purple", "black" = "red", "chinese" = "gold", "white" = "darkgreen")) +
    theme_minimal() +
    labs(
      title = paste(trait, "R² by Model and Ancestry"),
      x = "Model",
      y = "R²",
      fill = "Ancestry"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
})

# Save plots
plot_paths <- paste0("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/plots/barplot_preds_models_S_A_ancestries_", unique(summary_data_ancestries$Trait), ".png")
mapply(ggsave, filename = plot_paths, plot = plot_list, width = 8, height = 6)

# Display plots
print(plot_list)
