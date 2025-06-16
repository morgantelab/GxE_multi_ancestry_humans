rm(list = ls())

# Load necessary libraries
library(ggplot2)
library(readr)
library(data.table)
library(genio)
library(dplyr)

# Define plot save directory
plot_save_dir <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/plots/"

# --- Load Eigen Data ---
# PCrelate Eigenvectors
eigenvec_path_pcrelate <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/pca_for_pcrelate.rds"
eigenvec_data_pcrelate <- readRDS(eigenvec_path_pcrelate)

# GCTA Eigenvectors and Eigenvalues
eigenvec_path_gcta <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/pca_by_gcta.eigenvec"
eigenval_path_gcta <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/pca_by_gcta.eigenval"
eigenvec_data_gcta <- read_eigenvec(eigenvec_path_gcta)
eigenval_gcta <- read.table(eigenval_path_gcta)

# Struct MTD Eigenvectors
eigenvec_path_struct <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/eigen_by_struct.rds"
eigenvec_data_struct <- readRDS(eigenvec_path_struct)

# Load additional phenotype/ethnicity data
ethnicity_file_path <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/data2_20240930.RData"
load(ethnicity_file_path)
new_dataset <- dtt

# Merge with ethnicity data
# Extract rownames as IDs and convert to a data frame
eigenvec_pcrelate_df <- data.frame(ID = rownames(eigenvec_data_pcrelate$vectors),
                                   eigenvec_data_pcrelate$vectors, 
                                   stringsAsFactors = FALSE)
# Convert both IDs to character if necessary
eigenvec_pcrelate_df$ID <- as.character(eigenvec_pcrelate_df$ID)
new_dataset$ID <- as.character(new_dataset$ID)

# Merge the datasets by ID
merged_data_pcrelate <- merge(eigenvec_pcrelate_df, new_dataset, by = "ID")

# Extract IDs from the fam column
eigenvec_gcta_df <- data.frame(ID = eigenvec_data_gcta$fam$id, 
                               eigenvec_data_gcta$eigenvec, 
                               stringsAsFactors = FALSE)

# Convert GCTA IDs if necessary and merge
eigenvec_gcta_df$ID <- as.character(eigenvec_gcta_df$ID)
merged_data_gcta <- merge(eigenvec_gcta_df, new_dataset, by = "ID")

# Extract rownames as IDs and convert to a data frame
eigenvec_struct_df <- data.frame(ID = rownames(eigenvec_data_struct$vectors),
                                 eigenvec_data_struct$vectors, 
                                 stringsAsFactors = FALSE)

# Convert both IDs to character to ensure compatibility
eigenvec_struct_df$ID <- as.character(eigenvec_struct_df$ID)
merged_data_struct <- merge(eigenvec_struct_df, new_dataset, by = "ID")

colnames(merged_data_pcrelate)[2:ncol(eigenvec_pcrelate_df)] <- paste0("PC", 1:(ncol(eigenvec_pcrelate_df) - 1))
colnames(merged_data_gcta)[2:ncol(eigenvec_gcta_df)] <- paste0("PC", 1:(ncol(eigenvec_gcta_df) - 1))
colnames(merged_data_struct)[2:ncol(eigenvec_struct_df)] <- paste0("PC", 1:(ncol(eigenvec_struct_df) - 1))



# --- PC vs PC Plots ---
# Function to plot PC vs PC
plot_pcs <- function(data, pc_x, pc_y, title) {
  ggplot(data, aes_string(x = pc_x, y = pc_y, color = "ethn1_consolidated")) +
    geom_point() +
    theme_minimal() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(title = title, x = pc_x, y = pc_y, color = "Ethnicity")
}

# PCrelate PC vs PC plots
plot_pc1_pc2_pcrelate <- plot_pcs(merged_data_pcrelate, "PC1", "PC2", "PC1 vs PC2 (PCrelate)")
plot_pc1_pc3_pcrelate <- plot_pcs(merged_data_pcrelate, "PC1", "PC3", "PC1 vs PC3 (PCrelate)")
plot_pc2_pc3_pcrelate <- plot_pcs(merged_data_pcrelate, "PC2", "PC3", "PC2 vs PC3 (PCrelate)")

# GCTA PC vs PC plots
plot_pc1_pc2_gcta <- plot_pcs(merged_data_gcta, "PC1", "PC2", "PC1 vs PC2 (GCTA)")
plot_pc1_pc3_gcta <- plot_pcs(merged_data_gcta, "PC1", "PC3", "PC1 vs PC3 (GCTA)")
plot_pc2_pc3_gcta <- plot_pcs(merged_data_gcta, "PC2", "PC3", "PC2 vs PC3 (GCTA)")

# Struct MTD PC vs PC plots
plot_pc1_pc2_struct <- plot_pcs(merged_data_struct, "PC1", "PC2", "PC1 vs PC2 (Struct MTD)")
plot_pc1_pc3_struct <- plot_pcs(merged_data_struct, "PC1", "PC3", "PC1 vs PC3 (Struct MTD)")
plot_pc2_pc3_struct <- plot_pcs(merged_data_struct, "PC2", "PC3", "PC2 vs PC3 (Struct MTD)")

# --- Scree Plots ---
# PCrelate Scree plot
eigenval_pcrelate <- eigenvec_data_pcrelate$values
eigenval_df_pcrelate <- data.frame(PC = 1:length(eigenval_pcrelate), Variance = eigenval_pcrelate)
scree_plot_pcrelate <- ggplot(eigenval_df_pcrelate, aes(x = PC, y = Variance)) +
  geom_line() + geom_point() +
  theme_minimal() +
  labs(title = "Scree Plot (PCrelate)", x = "Principal Component", y = "Eigenvalue")

# GCTA Scree plot
eigenval_df_gcta <- data.frame(PC = 1:nrow(eigenval_gcta), Variance = eigenval_gcta$V1)
scree_plot_gcta <- ggplot(eigenval_df_gcta, aes(x = PC, y = Variance)) +
  geom_line() + geom_point() +
  theme_minimal() +
  labs(title = "Scree Plot (GCTA)", x = "Principal Component", y = "Eigenvalue")

# Struct MTD Scree plot
eigenval_struct <- eigenvec_data_struct$values
eigenval_df_struct <- data.frame(PC = 1:length(eigenval_struct), Variance = eigenval_struct)
scree_plot_struct <- ggplot(eigenval_df_struct, aes(x = PC, y = Variance)) +
  geom_line() + geom_point() +
  theme_minimal() +
  labs(title = "Scree Plot (Struct MTD)", x = "Principal Component", y = "Eigenvalue")

# --- Cumulative Scree Plots ---
# Function to create cumulative scree plot
cumulative_scree_plot <- function(eigenvalues, title) {
  total_variance <- sum(eigenvalues)
  cumulative_variance_explained <- cumsum(eigenvalues) / total_variance
  eigenval_df <- data.frame(PC = 1:length(cumulative_variance_explained), CumulativeVariance = cumulative_variance_explained)
  
  ggplot(eigenval_df, aes(x = PC, y = CumulativeVariance)) +
    geom_line() + geom_point() +
    theme_minimal() +
    labs(title = title, x = "Principal Component", y = "Cumulative Proportion of Variance Explained")
}

# PCrelate Cumulative Scree plot
cumulative_scree_plot_pcrelate <- cumulative_scree_plot(eigenval_pcrelate, "Cumulative Scree Plot (PCrelate)")

# GCTA Cumulative Scree plot
cumulative_scree_plot_gcta <- cumulative_scree_plot(eigenval_gcta$V1, "Cumulative Scree Plot (GCTA)")

# Struct MTD Cumulative Scree plot
cumulative_scree_plot_struct <- cumulative_scree_plot(eigenval_struct, "Cumulative Scree Plot (Struct MTD)")

# --- Save Plots to Files ---
# Save the PCrelate plots
ggsave(paste0(plot_save_dir, "pc1_vs_pc2_pcrelate.png"), plot_pc1_pc2_pcrelate, width = 8, height = 6)
ggsave(paste0(plot_save_dir, "pc1_vs_pc3_pcrelate.png"), plot_pc1_pc3_pcrelate, width = 8, height = 6)
ggsave(paste0(plot_save_dir, "pc2_vs_pc3_pcrelate.png"), plot_pc2_pc3_pcrelate, width = 8, height = 6)

# Save the GCTA plots
ggsave(paste0(plot_save_dir, "pc1_vs_pc2_gcta.png"), plot_pc1_pc2_gcta, width = 8, height = 6)
ggsave(paste0(plot_save_dir, "pc1_vs_pc3_gcta.png"), plot_pc1_pc3_gcta, width = 8, height = 6)
ggsave(paste0(plot_save_dir, "pc2_vs_pc3_gcta.png"), plot_pc2_pc3_gcta, width = 8, height = 6)

# Save the Struct MTD plots
ggsave(paste0(plot_save_dir, "pc1_vs_pc2_struct.png"), plot_pc1_pc2_struct, width = 8, height = 6)
ggsave(paste0(plot_save_dir, "pc1_vs_pc3_struct.png"), plot_pc1_pc3_struct, width = 8, height = 6)
ggsave(paste0(plot_save_dir, "pc2_vs_pc3_struct.png"), plot_pc2_pc3_struct, width = 8, height = 6)

# Save the scree plots
ggsave(paste0(plot_save_dir, "scree_plot_pcrelate.png"), scree_plot_pcrelate, width = 10, height = 6)
ggsave(paste0(plot_save_dir, "scree_plot_gcta.png"), scree_plot_gcta, width = 10, height = 6)
ggsave(paste0(plot_save_dir, "scree_plot_struct.png"), scree_plot_struct, width = 10, height = 6)

# Save cumulative scree plots
ggsave(paste0(plot_save_dir, "cumulative_scree_plot_pcrelate.png"), cumulative_scree_plot_pcrelate, width = 10, height = 6)
ggsave(paste0(plot_save_dir, "cumulative_scree_plot_gcta.png"), cumulative_scree_plot_gcta, width = 10, height = 6)
ggsave(paste0(plot_save_dir, "cumulative_scree_plot_struct.png"), cumulative_scree_plot_struct, width = 10, height = 6)

