# Load necessary libraries
library(ggplot2)
library(reshape2)
library(data.table)

# Set working directory
setwd("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/plots")

# Define lists of file paths for VCEm and varabs
VCEm_file_paths <- list(
  DP_plink = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/VCEm_DP_plink_G.csv",
  SP_plink = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/VCEm_SP_plink_G.csv",
  PP_plink = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/VCEm_PP_plink_G.csv",
  DP_pcrelate = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/VCEm_DP_pcrelate_G.csv",
  SP_pcrelate = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/VCEm_SP_pcrelate_G.csv",
  PP_pcrelate = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/VCEm_PP_pcrelate_G.csv"
)

varabs_file_paths <- list(
  DP_plink = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_DP_plink_G.csv",
  SP_plink = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_SP_plink_G.csv",
  PP_plink = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_PP_plink_G.csv",
  DP_pcrelate = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_DP_pcrelate_G.csv",
  SP_pcrelate = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_SP_pcrelate_G.csv",
  PP_pcrelate = "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/varabs_PP_pcrelate_G.csv"
)

# Define function to generate and save plots for each file
generate_and_save_plots <- function(VCEm_file_path, varabs_file_path, identifier) {
  
  # Read the VCEm data
  VCEm <- fread(VCEm_file_path)
  
  # Prepare VCEm data for plotting
  all_samples <- data.frame(
    Iteration = seq(1, nrow(VCEm)),
    Intercept = VCEm$int,
    GeneticVariance = VCEm$G,
    ResidualVariance = VCEm$res
  )
  
  # Melt the data for ggplot2
  VCEm_melted <- melt(all_samples, id.vars = "Iteration", variable.name = "Component", value.name = "Value")
  
  # Plot trace plots for intercept and variances
  intercept_variances_plot <- ggplot(VCEm_melted, aes(x = Iteration, y = Value, color = Component)) +
    geom_point(size = 0.5) +
    facet_wrap(~ Component, scales = "free_y") +
    ggtitle(paste("Trace Plots for Intercept and Variances -", identifier)) +
    xlab("Iteration") +
    ylab("Value")
  
  # Save the intercept and variances plot
  ggsave(paste0("G_trace_plots_intercept_and_variances_", identifier, ".png"), plot = intercept_variances_plot, width = 10, height = 6)
  
  # Read the varabs data
  varabs <- fread(varabs_file_path)
  
  # Prepare varabs data for plotting
  varabs_samples <- data.frame(
    Iteration = seq(1, nrow(varabs)),
    FixedEffectsVariance = varabs$V_X,
    AdditiveGeneticVariance = varabs$V_G
  )
  
  # Melt the varabs data for ggplot2
  varabs_melted <- melt(varabs_samples, id.vars = "Iteration", variable.name = "Component", value.name = "Variance")
  
  # Plot trace plots for the variances
  variances_plot <- ggplot(varabs_melted, aes(x = Iteration, y = Variance, color = Component)) +
    geom_point(size = 0.5) +
    facet_wrap(~ Component, scales = "free_y") +
    ggtitle(paste("Trace Plots for Variance Components -", identifier)) +
    xlab("Iteration") +
    ylab("Variance")
  
  # Save the variance components plot
  ggsave(paste0("G_trace_plots_variance_components_", identifier, ".png"), plot = variances_plot, width = 10, height = 6)
}

# Iterate over each pair of VCEm and varabs file paths
for (name in names(VCEm_file_paths)) {
  VCEm_file_path <- VCEm_file_paths[[name]]
  varabs_file_path <- varabs_file_paths[[name]]
  generate_and_save_plots(VCEm_file_path, varabs_file_path, name)
}
