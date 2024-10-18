# Load Required Libraries
library(tidyverse)
library(viridis)

setwd("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/plots")

load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/data2_20240930.RData")
dt <- dtt

# Define the ethnicity columns
ethnicity_columns <- c("ethn1_white", "ethn1_black", "ethn1_asian", "ethn1_mixed", "ethn1_chinese")

# List of environmental variables for which to create violin plots
env_variables <- c("Townsend", "act0_d", "TVtime", "sleep_d",
                   "veg_cook", "tea", "waist", "coffee", "BFP", 
                   "smoking_now", "fish_oily", "fish_lean", "meat_proc", "poultry",
                   "beef", "lamb", "pork", "cheese", "salt", "alc1", "getup", "smoked_past", "sleep_sd")

# Loop through each environmental variable
for (env_var in env_variables) {
  
  # Create a long dataframe for the current environmental variable across ethnicities and cohorts
  dt_long <- dt %>%
    pivot_longer(cols = all_of(ethnicity_columns), names_to = "Ethnicity", values_to = "Presence") %>%
    filter(Presence == 1) %>%
    select(coh, Ethnicity, all_of(env_var)) %>%
    rename(Value = all_of(env_var)) %>%
    mutate(
      Ethnicity = factor(Ethnicity, levels = ethnicity_columns, labels = c("White", "Black", "Asian", "Mixed", "Chinese")),
      Cohort = as.factor(coh)  # Convert 'coh' column to a factor
    )
  
  # Create the violin plot with box plot inside for the current environmental variable by cohort
  #  plot <- ggplot(dt_long, aes(x = Ethnicity, y = Value, fill = Ethnicity)) +
  #    geom_violin(trim = FALSE) +
  #    geom_boxplot(width = 0.1, fill = "white") +
  #    facet_wrap(~Cohort, scales = "free_y") +  # Facet by cohort
  #    scale_fill_viridis_d(begin = 0, end = 1, direction = 1) +
  #    labs(title = paste(env_var, "Distribution by Ethnicity and Cohort"), x = "", y = env_var) +
  #    theme_light() +
  #    theme(
  #      legend.position = "none",
  #      strip.text.x = element_text(size = 12, face = "bold"),
  #      axis.title.x = element_blank(),
  #      axis.title.y = element_blank(),
  #      axis.text.x = element_blank(),
  #      axis.ticks.x = element_blank()
  #    )
  
  # Create the violin plot with box plot inside for the current environmental variable by cohort
  plot <- ggplot(dt_long, aes(x = Ethnicity, y = Value, fill = Ethnicity)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.1, fill = "white") +
    scale_fill_viridis_d(begin = 0, end = 1, direction = 1) +
    facet_wrap(~Cohort, scales = "free_y", ncol = 4) +  # Facet by cohort
    labs(title = paste(env_var, "Distribution by Ethnicity and Cohort"), x = "", y = env_var) +
    theme_light() +
    theme(
      legend.position = "bottom",
      legend.justification = "right",
      legend.box.just = "right",
      legend.box.margin = margin(0, 0, 0, 0),
      strip.text.x = element_text(size = 12, face = "bold"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) 
  
  # Save the plot as a PNG file
  file_name <- paste0(env_var, "_ethnicity_cohort_violin_box_plot_improved_20241015.png")
  ggsave(file_name, plot, width = 8, height = 6, dpi = 300)
}

