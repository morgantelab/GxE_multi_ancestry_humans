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
folds <- 1:5  # Folds 1 to 5

# Generate the correct fold-based dataset with explicit Trait, Fold, and Model
fold_data <- expand.grid(Trait = traits, Fold = folds, Model = models, stringsAsFactors = FALSE) %>%
  mutate(File = sprintf("PREDs_%s_Fold_%d_%s.csv", Trait, Fold, Model))

# Merge with summary data to fill in R² values
summary_data_folds <- fold_data %>%
  left_join(summary_data, by = "File") %>%  # Merging to fill R² values
  select(Trait, Fold, Model, R_Squared)  # Keep only necessary columns

# Convert Fold to factor for proper ordering
summary_data_folds <- summary_data_folds %>%
  mutate(Fold = as.factor(Fold))

# Verify extraction correctness
print(table(summary_data_folds$Trait))  # Ensure Trait values are correct
print(table(summary_data_folds$Model))  # Ensure Model values are correct
print(table(summary_data_folds$Fold))   # Ensure Fold values are correct

# Define colors for each trait
trait_colors <- c("DP" = "limegreen", "SP" = "blue", "PP" = "orange")

# Define shapes for each fold
fold_shapes <- c(1, 2, 3, 4, 8)  # Shapes for Fold_1 to Fold_5

# Create boxplots with styling
plot_list <- lapply(unique(summary_data_folds$Trait), function(trait) {
  ggplot(summary_data_folds %>% filter(Trait == trait), aes(x = Model, y = R_Squared, fill = Trait)) +
    geom_boxplot(color = "black", alpha = 0.5) +  # Black outline, transparent fill
    geom_jitter(aes(shape = Fold), color = trait_colors[trait], size = 3, width = 0.2) +  # Jittered points
    scale_fill_manual(values = trait_colors) +
    scale_shape_manual(values = fold_shapes) +
    theme_minimal() +
    labs(
      title = paste(trait, "R² by Model"),
      x = "Model",
      y = "R²",
      fill = "Trait",
      shape = "Fold"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
})

# Save plots
plot_paths <- paste0("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/plots/boxplot_preds_models_S_A_folds_", unique(summary_data_folds$Trait), ".png")
mapply(ggsave, filename = plot_paths, plot = plot_list, width = 8, height = 6)

# Display plots
print(plot_list)
