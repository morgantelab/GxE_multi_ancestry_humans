library(ggplot2)
library(data.table)

# Load KS test results
ks_results <- fread("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/gwas_snp_env/ks_test_results_all_envs.txt")

# Convert P-values to -log10(P) for better visualization
ks_results[, logP := -log10(P_value)]

# Extract the trait from the file name (DP0s, SP0s, PP0s)
ks_results[, Trait := fifelse(grepl("DP0s", File), "DP0s",
                              fifelse(grepl("SP0s", File), "SP0s",
                                      fifelse(grepl("PP0s", File), "PP0s", NA_character_)))]

# Extract ancestry and fold labels from the file names
ks_results[, Ancestry := fifelse(grepl("asian", File), "Asian",
                                 fifelse(grepl("black", File), "Black",
                                         fifelse(grepl("white", File), "White",
                                                 fifelse(grepl("mixed", File), "Mixed",
                                                         fifelse(grepl("chinese", File), "Chinese", NA_character_)))))]

ks_results[, Fold := fifelse(grepl("_1\\.", File), "1",
                             fifelse(grepl("_2\\.", File), "2",
                                     fifelse(grepl("_3\\.", File), "3",
                                             fifelse(grepl("_4\\.", File), "4",
                                                     fifelse(grepl("_5\\.", File), "5", NA_character_)))))]

# Extract Environment (Interaction term) from the filename
ks_results[, Environment := gsub("gxe_|_[0-9]+.*", "", File)]  # Remove "gxe_" and fold suffix

# Function to generate improved Manhattan plots
plot_ks_manhattan <- function(data, category, category_name, trait) {
  # Ensure category is a factor for ordering
  data[[category]] <- factor(data[[category]], levels=unique(data[[category]]))
  
  # Manhattan-style plot with jitter, better label handling, and facetting
  plot <- ggplot(data, aes_string(x="Environment", y="D", color="logP")) +
    geom_point(size=3.5, alpha=0.8, position=position_jitter(width=0.2)) +  # Jitter points to avoid overlap
    scale_color_gradientn(colors=c("blue", "purple", "red"), limits=c(min(data$logP), max(data$logP))) +  # Better color contrast
    labs(title=paste("KS Test D-Statistic for", trait, "-", category_name),
         x="Gene-Environment Interaction",
         y="KS D-Statistic",
         color="-log10(P-value)") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle=35, hjust=1, vjust=1, size=8),  # Tilt X-axis labels
      panel.grid.major = element_blank(),  # Remove major grid
      panel.grid.minor = element_blank(),  # Remove minor grid
      legend.position="right",
      strip.text = element_text(size=10, face="bold")  # Improve facet label readability
    ) +
    facet_wrap(~Ancestry, scales="free_x")  # Facet by ancestry to separate data
  
  # Save the plot
  ggsave(filename=paste0("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/plots/manhattan_plot_KS_", 
                         category_name, "_", trait, ".png"), 
         plot=plot, width=12, height=6, dpi=300)
}

# Generate plots for each BP trait
for (trait in unique(ks_results$Trait)) {
  print(paste("Processing Trait:", trait))
  
  # Plot for Ancestries
  ancestry_data <- ks_results[!is.na(Ancestry) & Trait == trait]
  if (nrow(ancestry_data) > 0) {
    plot_ks_manhattan(ancestry_data, "Ancestry", "Ancestry Groups", trait)
  }
  
  # Plot for Cross-Validation Folds
  fold_data <- ks_results[!is.na(Fold) & Trait == trait]
  if (nrow(fold_data) > 0) {
    plot_ks_manhattan(fold_data, "Fold", "Cross-Validation Folds", trait)
  }
}

print("✅ Improved Manhattan plots for KS test statistics with environments saved!")



ggplot(ks_results, aes(x=D, y=logP, label=Interaction)) +
  geom_point(aes(color=logP), size=4, alpha=0.7) +
  geom_text(data=subset(ks_results, P_value < 1e-5), aes(label=Interaction), hjust=1.2, vjust=1.2, size=4) +
  scale_color_gradient(low="blue", high="red") +
  labs(title="KS Test: D-Statistic vs. -log10(P-value)",
       x="KS D-Statistic",
       y="-log10(P-value)",
       color="Significance") +
  theme_minimal()


library(reshape2)

# Convert data to matrix form for heatmap
ks_matrix <- dcast(ks_results, Interaction ~ File, value.var="D")

# Heatmap visualization
ggplot(melt(ks_matrix), aes(x=variable, y=Interaction, fill=value)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="red") +
  labs(title="Heatmap of KS D-Statistic",
       x="Result Files",
       y="Interaction Term",
       fill="D-Statistic") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1))

qqplot(qunif(ppoints(nrow(ks_results))), ks_results$D, main="QQ Plot of KS D-Statistic",
       xlab="Expected D-Statistic (Uniform)", ylab="Observed D-Statistic", pch=19, col="blue")
abline(0, 1, col="red", lty=2)  # Reference line

