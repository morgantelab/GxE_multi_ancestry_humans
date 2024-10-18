setwd("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/plots")

load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/data2_20240930.RData")


# Load Required Libraries
library(ggplot2)

# Function to create and save box plot for a given DP/SP column across ethnicities
create_boxplot <- function(data, value_column, ethnicity_columns, file_name) {
  # Pivot data to long format for ethnicity columns
  data_long <- data %>%
    pivot_longer(cols = all_of(ethnicity_columns), names_to = "Ethnicity", values_to = "Presence") %>%
    filter(Presence == 1) %>%
    select(Ethnicity, all_of(value_column)) %>%
    rename(Value = all_of(value_column)) %>%
    mutate(Ethnicity = factor(Ethnicity, levels = ethnicity_columns, labels = c("White", "Mixed", "Asian", "Black", "Chinese")))
  
  # Calculate the number of individuals per ethnicity
  n_counts <- data_long %>%
    group_by(Ethnicity) %>%
    summarize(n = n())
  
  # Create the box plot and add 'n' as text on the top-right of each box
  p <- ggplot(data_long, aes(x = Ethnicity, y = Value)) +
    geom_boxplot(fill = "blue", color = "black") +
    labs(title = paste(value_column, "Distribution Across Ethnicities"),
         x = "Ethnicity",
         y = value_column) +
    # Add number of individuals to the top of each box plot
    geom_text(data = n_counts, aes(x = Ethnicity, y = Inf, label = paste("n =", n)),
              vjust = -0.5, size = 5, color = "red")
  
  # Save the box plot
  ggsave(file_name, plot = p, width = 8, height = 6, dpi = 300)
}

# Assuming the data frame is named 'dtt' after loading the .RData file

# List of value columns (variables of interest)
value_columns <- c("DP0a", "SP0a", "PP0a")
# Ethnicity columns
ethnicity_columns <- c("ethn1_white", "ethn1_mixed", "ethn1_asian", "ethn1_black", "ethn1_chinese")

# Loop through each value column and create a single box plot across ethnicities
for (value_col in value_columns) {
  file_name <- paste("boxplot_of", value_col, "_across_ethnicities.png", sep = "")
  create_boxplot(dtt, value_col, ethnicity_columns, file_name)
}
