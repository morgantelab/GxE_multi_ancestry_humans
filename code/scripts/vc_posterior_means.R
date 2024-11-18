#clean enviroment
rm(list = ls())

library(dplyr)

# Define the list of files to process
file_paths <- list(
  DP = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_DP_plink_G_pcs.csv",
  SP = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_SP_plink_G_pcs.csv",
  PP = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_PP_plink_G_pcs.csv"
)

# Initialize a list to store the column means for each file
column_means_list <- list()

# Loop through each file
for (condition in names(file_paths)) {
  # Read the CSV file
  data <- read.csv(file_paths[[condition]])
  
  # Remove the column named "X" if it exists
  if ("X" %in% colnames(data)) {
    data <- data[ , !(colnames(data) %in% "X")]
  }
  
  # Calculate the mean of each column (excluding non-numeric columns)
  column_means <- colMeans(data, na.rm = TRUE)
  
  # Store the results in the list
  column_means_list[[condition]] <- column_means
  
  # Print the condition name and the column means
  cat("\nCondition:", condition, "\n")
  print(column_means)
}

# Combine results into a single data frame
#combined_means <- do.call(rbind, column_means_list)
#combined_means <- as.data.frame(combined_means)
#combined_means$Condition <- rownames(combined_means)
#rownames(combined_means) <- NULL

# Reorder the columns to place "Condition" as the first column
#combined_means <- combined_means[, c("Condition", setdiff(names(combined_means), "Condition"))]

# Define the path for the output CSV file
#output_file <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/combined_column_means.csv"

# Write the combined means to a CSV file
#write.csv(combined_means, file = output_file, row.names = FALSE)

# Print a message indicating the file was saved
#cat("Combined column means saved to:", output_file, "\n")

#combined_column_means <- combined_column_means %>% mutate(Condition = ifelse(Condition == "PP", "PP_Plink_X_separated", Condition))
