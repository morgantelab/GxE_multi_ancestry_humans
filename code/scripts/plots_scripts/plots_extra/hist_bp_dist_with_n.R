setwd("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/plots")

load("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/data2_20240930.RData")

# Load Required Libraries
library(ggplot2)

# Function to create and save histogram for a given DP/SP column and ethnicity
create_histogram <- function(data, value_column, ethnicity_column, file_name) {
  # Filter data for the ethnicity
  data_filtered <- data[data[[ethnicity_column]] == 1, ]
  
  # Calculate the number of individuals (rows)
  n <- nrow(data_filtered)
  
  # Create the histogram and add 'n' to the top-right
  p <- ggplot(data_filtered, aes_string(x = value_column)) +
    geom_histogram(bins = 30, fill = "blue", color = "black") +
    labs(title = paste(value_column, "for", sub("ethn1_", "", ethnicity_column)), # Shorten title
         x = value_column,
         y = "Frequency") +
    annotate("text", x = Inf, y = Inf, label = paste("n =", n), 
             hjust = 1.1, vjust = 1.5, size = 5, color = "red")  # Adjust the text position
  
  # Save the histogram
  ggsave(file_name, plot = p, width = 8, height = 6, dpi = 300)
}

# Assuming the data frame is named 'dtt' after loading the .RData file

# List of value columns and ethnicity columns
value_columns <- c("DP0a", "SP0a", "PP0a")
ethnicity_columns <- c("ethn1_white", "ethn1_mixed", "ethn1_asian", "ethn1_black", "ethn1_chinese")

# Loop through each combination of value and ethnicity columns
for (value_col in value_columns) {
  for (ethn_col in ethnicity_columns) {
    file_name <- paste("histogram_of", value_col, "vs", ethn_col, ".png", sep = "")
    create_histogram(dtt, value_col, ethn_col, file_name)
  }
}
