# all figures for manuscript


# Figure 1: raw data + counts ----
if(1){
  print('Fig 1')
  set.seed(1123)

  # Load necessary libraries
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(gridExtra)
  library(grid)

  # Load dataset
  file_path <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/subset_dataset_figure1.csv"
  df <- read_csv(file_path)

  # Capitalize only the first letter of each ethnicity value
  df$ethn1_consolidated <- paste0(toupper(substring(df$ethn1_consolidated, 1, 1)),
                                  tolower(substring(df$ethn1_consolidated, 2)))

  # Define ethnicity color mapping
  ethnicity_colors <- c(
    "White" = "limegreen",
    "Black" = "orange",
    "Asian" = "yellow",
    "Mixed" = "pink",
    "Chinese" = "red"
  )

  # Define common theme adjustments (with plot box border)
  common_theme <- theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.8),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.margin = margin(10, 10, 10, 10)  # Added margin to prevent cropping
    )

  # Define the order of ethnicity categories
  ethnicity_order <- c("Asian", "White", "Black", "Mixed", "Chinese")

  # Create figure (A) - Diastolic Blood Pressure (DP0a)
  fig_a <- ggplot(df, aes(x = factor(ethn1_consolidated, levels = ethnicity_order), y = DP0a, fill = ethn1_consolidated)) +
    geom_boxplot(color = "black", alpha = 0.7) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3.5,
                 fill = "deepskyblue3", color = "black") +
    scale_fill_manual(values = ethnicity_colors) +
    scale_x_discrete(labels = ethnicity_order) +
    common_theme +
    labs(y = "mmHg", x = "")
  # , tag='A')+
  #   theme(plot.tag = element_text(size=tagsize))

  # Create figure (B) - Systolic Blood Pressure (SP0a)
  fig_b <- ggplot(df, aes(x = factor(ethn1_consolidated, levels = ethnicity_order), y = SP0a, fill = ethn1_consolidated)) +
    geom_boxplot(color = "black", alpha = 0.7) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3.5,
                 fill = "deepskyblue3", color = "black") +
    scale_fill_manual(values = ethnicity_colors) +
    scale_x_discrete(labels = ethnicity_order) +
    common_theme +
    labs(y = "mmHg", x = "")
  # , tag='B')+
  #   theme(plot.tag = element_text(size=tagsize))

  # Create figure (C) - Pulse Pressure (PP0a) with improved y-axis
  fig_c <- ggplot(df, aes(x = factor(ethn1_consolidated, levels = ethnicity_order), y = PP0a, fill = ethn1_consolidated)) +
    geom_boxplot(color = "black", alpha = 0.7) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3.5,
                 fill = "deepskyblue3", color = "black") +
    scale_fill_manual(values = ethnicity_colors) +
    scale_x_discrete(labels = ethnicity_order) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +  # improved spacing
    common_theme +
    labs(y = "mmHg", x = "")
  # , tag='C')+
  #   theme(plot.tag = element_text(size=tagsize))

  # Create figure (D) - Ethnicity Count Bar Plot (with flexible y-axis)
  ethnicity_counts <- df %>%
    group_by(ethn1_consolidated) %>%
    summarise(count = n())

  fig_d <- ggplot(ethnicity_counts, aes(x = factor(ethn1_consolidated, levels = ethnicity_order), y = count, fill = ethn1_consolidated)) +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = ethnicity_colors) +
    scale_x_discrete(labels = ethnicity_order) +
    common_theme +
    labs(y = "Count", x = "")
  # , tag='D')+
  #   theme(plot.tag = element_text(size=tagsize))

  # Function to add labels (A, B, C, D) **INSIDE** the plot
  add_label <- function(plot, label) {
    arrangeGrob(
      textGrob(label, x = 0.1, y = 0.50, just = "left", gp = gpar(fontsize = 14, fontface = "bold")),
      plot,
      heights = c(0.05, 1),  # Reduced label height
      ncol = 1
    )
  }

  final_plot <- plot_grid(fig_a, fig_b, fig_c, fig_d, labels = c("A", "B", 'C', 'D'), ncol=2)
  final_plot
  root <- '/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/output/figs/final/'
  fname <- 'Figure_1.pdf'
  onepath <- paste0(root, fname)

  # Save final plot
  ggsave(plot = final_plot,
         filename = onepath,
         width = 17.4, height = 17.4, units = 'cm', dpi = 300, device = 'pdf')


  # Arrange figures with properly positioned labels **INSIDE IMAGE BOUNDS**
  # final_plot <- grid.arrange(
  #   add_label(fig_a, "A"),
  #   add_label(fig_b, "B"),
  #   add_label(fig_c, "C"),
  #   add_label(fig_d, "D"),
  #   ncol = 2, nrow = 2
  # )

  # # Save all figures in one PDF with **STRICT SIZE**
  # pdf ("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/plots/Figure_1.pdf", width = 8, height = 8)
  # grid.arrange(
  #   add_label(fig_a, "A"),
  #   add_label(fig_b, "B"),
  #   add_label(fig_c, "C"),
  #   add_label(fig_d, "D"),
  #   ncol = 2, nrow = 2
  # )
  # dev.off()

}

# Figure 2: PCs by method, ERM PCs ----
if(1){
  library(data.table)
  print('Fig 2')
  base_theme <- theme_bw(base_size = 12) +  # Bigger overall text
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(color = "black", size = 8),
      axis.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      plot.title = element_text(size = 8)
    )


  # data prep for MASSIVE eigen vector 24k x 24k files
  if(0){
    # Load lifestyle PCA Data
    eigenvec_data_lifestyle <- readRDS("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/E_eigen.rds")


    # Merge with ethnicity data
    eigenvec_lifestyle_df <- data.table(ID = rownames(eigenvec_data_lifestyle$vectors),
                                        eigenvec_data_lifestyle$vectors)

    # saveRDS(eigenvec_lifestyle_df, '/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/eigen_lifestyle_reformat.rds')
    # names( eigenvec_lifestyle_df)[1:4]

    # Load phenotype/ethnicity data
    rawdtt <- readRDS("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/data3_20250514.rds")

    # Ensure IDs are characters
    eigenvec_lifestyle_df$ID <- as.character(eigenvec_lifestyle_df$ID)

    dtt <- rawdtt[,c('ID', 'ethn1_consolidated')]
    dtt$ID <- as.character(dtt$ID)

    # Merge datasets by ID
    merged_data_lifestyle <- merge(eigenvec_lifestyle_df, dtt, by = "ID")

    # Rename PCs clearly
    colnames(merged_data_lifestyle)[2:ncol(eigenvec_lifestyle_df)] <- paste0("PC", 1:(ncol(eigenvec_lifestyle_df) - 1))
    # names(merged_data_lifestyle)[1:4]
    # length(merged_data_lifestyle)
    # names(merged_data_lifestyle)[24478:24562]
    #

    colkeep <- c('ID', 'PC1', 'PC2', 'ethn1_consolidated')
    final_merged_lifestyle_data <- data.table(merged_data_lifestyle[,colkeep,with=0])
    saveRDS(final_merged_lifestyle_data,  '/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/eigen_lifestyle_reformat.rds')






    # Panel A: PLINK PCA
    eigenvec_data_plink <- readRDS("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/pca_for_plink.rds")
    eigenvec_plink_df <- data.frame(ID = rownames(eigenvec_data_plink$vectors),
                                    eigenvec_data_plink$vectors,
                                    stringsAsFactors = FALSE)
    eigenvec_plink_df$ID <- as.character(eigenvec_plink_df$ID)
    merged_data_plink <- merge(eigenvec_plink_df, dtt, by = "ID")
    colnames(merged_data_plink)[2:ncol(eigenvec_plink_df)] <- paste0("PC", 1:(ncol(eigenvec_plink_df) - 1))

    colkeep <- c('ID', 'PC1', 'PC2', 'ethn1_consolidated')
    final_merged_plink_data <- data.table(merged_data_plink[,colkeep])
    saveRDS(final_merged_plink_data,  '/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/eigen_plink_reformat.rds')

    # Panel B: PC-Relate PCA
    eigenvec_data_pcrelate <- readRDS("/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/filtered_chr/pca_for_pcrelate.rds")
    eigenvec_pcrelate_df <- data.frame(ID = rownames(eigenvec_data_pcrelate$vectors),
                                       eigenvec_data_pcrelate$vectors,
                                       stringsAsFactors = FALSE)
    eigenvec_pcrelate_df$ID <- as.character(eigenvec_pcrelate_df$ID)
    merged_data_pcrelate <- merge(eigenvec_pcrelate_df, dtt, by = "ID")
    colnames(merged_data_pcrelate)[2:ncol(eigenvec_pcrelate_df)] <- paste0("PC", 1:(ncol(eigenvec_pcrelate_df) - 1))

    colkeep <- c('ID', 'PC1', 'PC2', 'ethn1_consolidated')
    final_merged_pcrelate_data <- data.table(merged_data_pcrelate[,colkeep])
    saveRDS(final_merged_pcrelate_data,  '/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/eigen_pcrelate_reformat.rds')


  }

  merged_data_pcrelate <- readRDS('data/eigen_pcrelate_reformat.rds')
  merged_data_plink <- readRDS('data/eigen_plink_reformat.rds')
  merged_data_lifestyle <- readRDS('data/eigen_lifestyle_reformat.rds')


  ethnicity_colors <- c(
    "White" = "limegreen",
    "Black" = "orange",
    "Asian" = "yellow",
    "Mixed" = "pink",
    "Chinese" = "red")

  # Plot panels
  mysize <- 1
  myalpha <- 0.7
  plot_plink <- ggplot(merged_data_plink, aes(x = PC1, y = PC2, color = ethn1_consolidated)) +
    geom_point(alpha = myalpha, size = mysize) +
    base_theme +
    theme(legend.position = "none") +
    scale_color_manual(values = ethnicity_colors) +
    labs(x = "Principal Component 1", y = "Principal Component 2")+
    labs(color='Ethnic Group')

  plot_pcrelate_no_legend <- ggplot(merged_data_pcrelate, aes(x = PC1, y = PC2, color = ethn1_consolidated)) +
    geom_point(alpha = myalpha, size = mysize) +
    base_theme +
    theme(legend.position = "none") +
    scale_color_manual(values = ethnicity_colors) +
    labs(x = "Principal Component 1", y = "Principal Component 2")

  pc1_vs_pc2_plot <- ggplot(merged_data_lifestyle, aes(x = PC1, y = PC2, color = ethn1_consolidated)) +
    geom_point(alpha = myalpha, size = mysize) +
    base_theme +
    theme(legend.position = "none") +
    scale_color_manual(values = ethnicity_colors) +
    labs(x = "Principal Component 1",
         y = "Principal Component 2")

  # Scraped legend
  nline <- 1.2
  plot_for_legend <- ggplot(merged_data_pcrelate, aes(x = PC1, y = PC2, color = ethn1_consolidated)) +
    geom_point(size = 3) +
    theme_minimal(base_size = 11) +
    guides(color=guide_legend("Ethnicity"))+
    theme(
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 11),
      legend.key.height = unit(nline, "lines")  # ⬅️ more vertical spacing
    ) +
    scale_color_manual(
      values = ethnicity_colors,
      guide = guide_legend(
        override.aes = list(size = 1.5),  # ⬅️ larger points in legend
        keyheight = unit(nline, "lines")  # ⬅️ match vertical space
      )
    )

  legend_shared <- get_legend(plot_for_legend)


  final_plot <- plot_grid(plot_plink,
                          plot_pcrelate_no_legend,
                          pc1_vs_pc2_plot,
                          legend_shared,
                          labels = c("A", "B", 'C', ' '),
                          ncol = 2)

  final_plot

  root <- '/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/output/figs/final/'
  fname <- 'Figure_2.pdf'
  twopath <- paste0(root, fname)

  # Save final plot
  ggsave(plot = final_plot,
         filename = twopath,
         width = 17.4, height = 17.4, units = 'cm', dpi = 300, device = 'pdf')
}

# Figure 3: Variance component comparison ----
if(1){
  print('Fig 3')
  library(data.table)
  library(ggplot2)
  library(ggpubr)


  # full table of variance component effect sizes
  # matte_old <- fread('output/vault/fullport/varcomps.txt')
  matte <- fread('data/nklimko/new_varcomp_1-7.txt')

  matte$comp2 <- matte$comp

  matte$comp2 <- gsub('^X1$','A_X1',matte$comp2)
  matte$comp2 <- gsub('^X2$','B_X2',matte$comp2)
  matte$comp2 <- gsub('^E$','C_E',matte$comp2)
  matte$comp2 <- gsub('^G$','D_G',matte$comp2)
  matte$comp2 <- gsub('^m7$','E_M7',matte$comp2)
  matte$comp2 <- gsub('^GxE$','F_GxE',matte$comp2)

  # matte$fax <- paste0(matte$comp,'_from_',matte$model)
  matte$fax <- paste0(matte$model, '_', matte$comp2)
  matte$modnum <- gsub('Model ', '', matte$model)

  matte <- data.table(matte)

  mat_dp <- matte[trait=='DP',]
  mat_sp <- matte[trait=='SP',]
  mat_pp <- matte[trait=='PP',]

  dp_bind <- mat_dp[,sum(mean), by=.(model)]
  sp_bind <- mat_sp[,sum(mean), by=.(model)]
  pp_bind <- mat_pp[,sum(mean), by=.(model)]

  remolin <- function(datanow, fulltitle, montag=NULL, mod7=TRUE, size_point=4, mjust=0, textjust=0, buffer=0, custard = 'nerg'){
    plothole <- ggplot(data=datanow, aes(x=model, y=mean))+
      scale_color_manual(labels = c('Demographics', 'Structure', 'Lifestyle', 'Genetics', 'Ancestry-by-lifestyle', 'Genetics-by-lifestyle', ' '), values= c("#E69F00","#56B4E9","#009E73","purple",'pink', "brown",'black'))+
      geom_point(aes(color=comp2), position = position_dodge2(width = 1, preserve = "single"), size = size_point/2) +
      geom_errorbar(aes(ymin = lower, ymax = upper, color=comp2),
                    width = 1, linewidth = 0.4, alpha = 1, position = position_dodge2(width = 1, preserve = "single")) +
      geom_text(data=custard, aes(x = model, y = textjust, label= round(V1, 2)), color='black', show.legend = FALSE)+
      geom_text(aes(x = 'Model 2', y = textjust + buffer, label=" ", color='white', fill='white'), show.legend = FALSE)+
      labs(title=fulltitle, color=NULL, tag=montag)+
      xlab(NULL)+
      ylab('Variance')+
      theme_minimal()+
      theme(legend.position = 'bottom',
            # axis.text.x = element_text(hjust = mjust),
            # axis.text.x = element_text(size=8),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))+
      # geom_vline(xintercept = verty, linetype=2, color='lightgray')+
      # geom_hline(yintercept = 0, color='black')+
      # geom_text(data=modcode, hjust='middle', aes(label=mn, x = x, y = y), inherit.aes=FALSE)+
      guides(color = guide_legend(override.aes = list(size = size_point / 2)))


    return(plothole)
  }

  dpom <- remolin(mat_dp, fulltitle=NULL, mod7=TRUE, textjust=23, buffer=2, custard=dp_bind)
  spom <- remolin(mat_sp, fulltitle=NULL, mod7=TRUE, textjust=22, buffer=2, custard=sp_bind)
  ppom <- remolin(mat_pp, fulltitle=NULL, mod7=TRUE, textjust=25, buffer=2, custard=pp_bind)

  final_plot <- ggarrange(dpom, spom, ppom, ncol=1, labels = c("A", "B", 'C'), common.legend = TRUE, legend='bottom')
  final_plot

  root <- '/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/output/figs/final/'
  fname <- 'Figure_3.pdf'
  threepath <- paste0(root, fname)

  # Save final plot
  ggsave(plot = final_plot,
         filename = threepath,
         width = 17.4, height = 17.4, units = 'cm', dpi = 300, device = 'pdf')


}

# Figure 4: cross-ethnic prediction accuracy (bar plots x3) ----
if(1){
  print('Fig 4')
  # rm(list=ls()); gc()
  set.seed(1123)

  # Load required libraries
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(cowplot)
  library(scales)
  library(ggpubr)

  size4 <- 10

  # Load the summary file
  summary_file <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/all_model_summary_ethn_corrected_use.csv"
  summary_data <- read_csv(summary_file)

  # summary_data <- data.table(summary_data)
  # Define levels
  traits <- c("SP", "DP", "PP")
  ancestries <- c("asian", "mixed", "black", "chinese", "white")
  model_levels <- c("X1", "X1_X2", "X1_X2_G", "X1_X2_E", "X1_X2_G_E", "X1_X2_G_E_GE")
  model_labels <- c("Model 1", "Model 2", "Model 3", "Model 4", "Model 5", "Model 6")

  # Filter and format data
  summary_data_ancestries <- summary_data %>%
    dplyr::filter(Trait %in% traits, Model %in% model_levels)  %>%
    mutate(
      Trait = factor(Trait, levels = traits),
      Ancestry = factor(str_to_title(Ancestry), levels = str_to_title(ancestries)),
      Model = factor(Model, levels = model_levels),
      models = factor(models, levels = model_labels)
    )

  # Define custom model colors to match the image
  model_colors <- c(
    "X1" = "#1b9e77",        # Teal/Dark Green (Model 1)
    "X1_X2" = "#d95f02",       # Orange (Model 2)
    "X1_X2_G" = "#7570b3",    # Purple (Model 3)
    "X1_X2_E" = "#e7298a",    # Pink/Magenta (Model 4)
    "X1_X2_G_E" = "#66a61e",   # Lighter Green (Model 5) - adjusted to be distinct from M1, based on the image
    "X1_X2_G_E_GE" = "#e6ab02" # Yellow/Gold (Model 6)
  )

  # Trait plot labels
  traits_titles <- c("DP" = "A", "SP" = "B", "PP" = "C")

  # Create individual bar plots with flexible y-axis
  plot_list <- lapply(names(traits_titles), function(trait) {
    ggplot(summary_data_ancestries %>% dplyr::filter(Trait == trait),
           aes(x = Ancestry, y = R_squared, fill = Model)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.8),
               color = "black", width = 0.7) +
      scale_fill_manual(
        values = model_colors,
        breaks = model_levels,
        labels = model_labels
      ) +
      scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
      theme_bw(base_size = size4 / 2) +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(color = "black", size = size4), # Increased size of x-axis labels
        axis.text.y = element_text(color = "black", size = size4), # Increased size of y-axis labels
        axis.title.y = element_text(face = "bold",  size = size4), # Increased size of y-axis title to match labels
        axis.title.x = element_blank(),
        legend.position = "right",
        # legend.position = "right",
        # legend.title = element_blank(),
        legend.title = element_text(size=size4),
        legend.text = element_text(size = size4), # Increased size of legend text
        legend.key.size = unit(0.4, "cm"),
        plot.title.position = "plot",
        plot.title = element_text(hjust = 0, face = "bold", size = size4 * 1.5, color = "black"), # Adjusted hjust to 0 for more left alignment
        panel.border = element_rect(colour = "black", fill=NA, size=0.5) # Increased size of box outline
      ) +
      labs(
        # title = traits_titles[[trait]],
        y = expression(italic(R^2)),
        fill='Model'
      )
  })

  # Stack the plots vertically
  final_plot <- ggarrange(plotlist = plot_list, ncol=1, labels = c("A", "B", 'C'), legend = 'right', common.legend = TRUE)

  final_plot

  root <- '/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/output/figs/final/'
  fname <- 'Figure_4.pdf'
  fourpath <- paste0(root, fname)

  # Save the final PDF
  ggsave(final_plot, filename = fourpath, height = 17.4, width = 17.4, units='cm', dpi = 300, device = "pdf")

}


# Figure 5: Kfold CV prediction accuracy (whiskers m1-6) ----
if(1){
  print('Fig 5')
  # Clear environment
  # rm(list = ls()); gc()
  set.seed(1123)

  # Load required libraries
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(cowplot)

  # Load updated data
  summary_file <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/data/model/all_model_summary_folds.csv"
  summary_data <- read_csv(summary_file)

  # Filter out GEselect rows
  summary_data <- summary_data %>% dplyr::filter(!grepl("GEselect", Model))

  # Define parameters
  traits <- c("SP", "DP", "PP")
  models <- c("X1", "X1_X2", "X1_X2_G", "X1_X2_E", "X1_X2_G_E", "X1_X2_G_E_GE")
  model_labels <- c("X1" = "Model 1", "X1_X2" = "Model 2", "X1_X2_G" = "Model 3",
                    "X1_X2_E" = "Model 4", "X1_X2_G_E" = "Model 5", "X1_X2_G_E_GE" = "Model 6")
  folds <- 1:5

  # Format data
  summary_data_folds <- summary_data %>%
    dplyr::filter(Trait %in% traits, Model %in% models) %>%
    mutate(
      Trait = factor(Trait, levels = traits),
      Fold = factor(Fold, levels = folds),
      Model = factor(Model, levels = models)
    )


  # Colorblind-friendly color palette (Set2)
  trait_colors <- c("DP" = "#66C2A5", "SP" = "#8DA0CB", "PP" = "#FC8D62")
  traits_titles <- c("DP" = "A", "SP" = "B", "PP" = "C")


  size5 <- 10

  # Generate boxplots with flexible y-axis per trait
  plot_list <- lapply(names(traits_titles), function(trait) {
    data_trait <- summary_data_folds %>% dplyr::filter(Trait == trait)

    # Calculate mean R² values for diamonds
    mean_data <- data_trait %>%
      group_by(Model) %>%
      summarize(mean_r2 = mean(R_squared, na.rm = TRUE), .groups = 'drop')

    p <- ggplot(data_trait, aes(x = Model, y = R_squared)) +
      geom_boxplot(color = "black", alpha = 0.6, fill = trait_colors[trait], width = 0.6, fatten=1) +
      geom_point(data = mean_data, aes(x = Model, y = mean_r2),
                 shape = 23, size = size5 / 5, fill = "#FFD700", color = "black", stroke = 0.6) +  # Filled golden diamonds # size 4.5
      scale_x_discrete(labels = model_labels) +
      theme_bw(base_size = size5) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_blank(),
        axis.text = element_text(color = "black", size = size5), # 28
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = size5), # 30
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, margin = margin(t = 12)),
        legend.position = "none",
        plot.title = element_text(hjust = -0.1, face = "bold")
      ) +
      labs(
        # title = traits_titles[[trait]],
        y = expression(italic(R^2))
      )

    # Custom y-axis for PP and SP
    if (trait == "PP") {
      p <- p + scale_y_continuous(
        limits = c(0.19, 0.235),
        breaks = seq(0.19, 0.235, by = 0.01),
        labels = scales::number_format(accuracy = 0.01)
      )
    } else if (trait == "SP") {
      # Auto-determine useful breaks with gap of 0.03
      y_min <- floor(min(data_trait$R_squared, na.rm = TRUE) * 50) / 50  # closest lower multiple of 0.02
      y_max <- ceiling(max(data_trait$R_squared, na.rm = TRUE) * 50) / 50  # closest upper multiple of 0.02
      p <- p + scale_y_continuous(
        breaks = seq(y_min, y_max, by = 0.03),
        labels = scales::number_format(accuracy = 0.01)
      )
    } else {
      p <- p + scale_y_continuous(labels = scales::number_format(accuracy = 0.01))
    }

    return(p)
  })

  # Combine plots into a column
  final_plot <- plot_grid(plotlist = plot_list, labels = c("A", "B", 'C'), ncol = 1)

  final_plot


  root <- '/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/output/figs/final/'
  fname <- 'Figure_5.pdf'
  fivepath <- paste0(root, fname)

  # Save final plot
  ggsave(plot = final_plot,
         filename = fivepath,
         width = 17.4, height = 17.4, units = 'cm', dpi = 300, device = 'pdf')
}



# Fig S1: Pop count ----
if(1){

  library(ggplot2)
  library(ggpubr)

  invr3 <- fread('output/vault/sups/invr3.tx')
  r2 <- fread('output/vault/sups/waterfall.txt')
  tr2 <- fread('output/vault/sups/waterfall_transpose.txt')

  molten <- melt(invr3)

  colorder <- names(r2)[41:2]


  key <- data.table(mikey=101:140, variable=colorder)
  key$mikey <- as.factor(key$mikey)

  names(key) <- c('mikey', 'Step')

  lava <- molten[key, on=.(Step)]

  lava2 <- lava[variable!='White' & variable!='Total' & variable!='Other/NA',]
  lava3 <- lava[variable=='White',]


  biglabel<- c('101' = 'base',
               '102' = 'genotyped',
               '103' = 'DP_SP',
               '104' = 'DP_SP_4sd',
               '105' = 'med_history',
               '106' = 'self_sex',
               '107' = 'age_pheno',
               '108' = 'age_recruit',
               '109' = 'townsend',
               '110' = 'act0_d',
               '111' = 'TVtime',
               '112' = 'sleep',
               '113' = 'smoke_current',
               '114' = 'veg_cook',
               '115' = 'fish_oily',
               '116' = 'fish_lean',
               '117' = 'meat_proc',
               '118' = 'poultry',
               '119' = 'beef',
               '120' = 'lamb',
               '121' = 'pork',
               '122' = 'cheese',
               '123' = 'salt',
               '124' = 'tea',
               '125' = 'alcohol',
               '126' = 'waist',
               '127' = 'mediblood_again',
               '128' = 'getup',
               '129' = 'coffee',
               '130' = 'smoke_past',
               '131' = 'bfp',
               '132' = 'age_bound',
               '133' = 'ethn_na_removed',
               '134' = 'ethn_other_removed',
               '135' = 'coffee_99',
               '136' = 'tea_99',
               '137' = 'veg_cook_99',
               '138' = 'tvtime_99',
               '139' = 'sleep_01',
               '140' = 'sleep_99')


  biglabel <- gsub('_',' ',biglabel)

  nubs <- data.table(rep(names(tr2), 2), c(rep(3, 7), rep(38, 7)), c(502460,9880,8060,1573,2956,472656,7335,
                                                                     321388,4957,3979,774,1746,309932,0))
  names(nubs) <- c('ethn', 'colm', 'count')
  nubs2 <- nubs[ethn!='White' & ethn!='variable' & ethn!='zz_other',]


  # nubs3$colm <- c(4,37)


  nubs3 <- nubs[ethn=='White',]

  nubs3$colm <- c(4,37)

  nubs3$countspace <- c(475000,305000)

  ethnicity_colors <- c(
    "White" = "limegreen",
    "Black" = "orange",
    "Asian" = "yellow",
    "Mixed" = "pink",
    "Chinese" = "red"
  )

  sqar <- 17.4
  sf1 <- 5.9

  linesize <- 0.8
  waterlayers <- list(scale_color_manual(values=ethnicity_colors, drop=FALSE),
                      theme_minimal(),
                      geom_step(aes(group=variable), size=linesize, show.legend=TRUE),
                      scale_x_discrete("Stage",labels=biglabel),
                      theme(axis.text.x = element_text(angle=90),
                            legend.text = element_text(size=sf1 - 1),
                            legend.title = element_text(size=sf1),
                            axis.text = element_text(size=sf1),
                            axis.title = element_text(size=sf1 + 1),
                            legend.position = 'bottom',
                            panel.grid.major.x=element_blank()),
                      labs(y = 'Individuals', color='Ethnic Group')
  )


  plok1 <- ggplot(lava2, aes(x=mikey, y = value, color=variable))+
    waterlayers+
    geom_text(data=nubs2, size = sf1 / 2, aes(x=colm, y=count + 150, label=count), inherit.aes=FALSE)+
    labs(tag='B')

  plok2 <- ggplot(lava3, aes(x=mikey, y = value, color=variable))+
    waterlayers+
    geom_text(data=nubs3, size = sf1 / 2, aes(x=colm, y=countspace, label=count), inherit.aes=FALSE)+
    labs(tag='A')

  # plokky <- plot_grid(plok2, plok1, labels =c('A', 'B'), ncol=2)

  plokky <- ggarrange(plok2, plok1, ncol = 2, common.legend = TRUE, legend = 'bottom')

  ggsave(plokky, height=sqar / 1.5, width=sqar, units='cm', device='png', filename='output/figs/final/Figure_S1_alt.png', dpi=300)
}

# Fig S2: Env covar ----
if(1){

  library(ggcorrplot)
  load('data/Emat_20250514.RData')
  # dim(Emat)
  d2 <- cov(Emat)
  golden <- ggcorrplot(d2, type='upper', hc.order=TRUE, show.diag = TRUE)

  # golden

  sqar <- 17.4
  ggsave(golden, height=sqar, width=sqar, units='cm', device='png', filename='output/figs/final/Figure_S2.png', dpi=300)


}


# S3-S5: QQ Plots ----
if(1){

  # needed for toTitleCase in title remake
  library(tools)
  library(data.table)
  library(ggplot2)
  library(ggpubr)


  claude <- function(trait, env){

    #enforce upper case
    trait <- toupper(trait)

    #read in data
    newp <- readRDS(paste0('output/gpcxe/tables_gxegwas/', trait,'/',env,'.rds'))

    #convert to numeric
    newp$P <- as.numeric(newp$P)

    # Remake title
    # newtitle <- paste0(trait, ' ' ,env, ' new')
    env2 <- gsub('_',' ',env)
    newtitle <- toTitleCase(env2)

    newcut <- 0
    minder <- paste0(trait,'_',env)
    # print(minder)

    if(nrow(mindful[merged==minder,]) > 0) {
      # if (minder %in% mindful$cutoff){
      newcut <- mindful[merged==minder, cutoff]
    }
    # print(newcut)

    #qqlim funciton call
    newhead <- qqlim(unlist(newp$P), fulltitle=newtitle, saving=0, size_point=0.75, fullcut=0, cutoff=newcut)

    return(newhead)

  }

  qqlim <- function(rawvector, outpath=NA, file=NA, size_point=0.75, size_lamb=4, saving=0, fulltitle=NULL, fullcut=0, cutoff=0){

    # data filtering
    pvector <- na.omit(rawvector)

    # Calculate lambda
    chisq <- qchisq(1-pvector, 1)
    GIF <- median(chisq)/qchisq(0.5,1)
    smallGIF <- round(GIF, 3)

    # expected and observed
    expected <- -log10(ppoints(length(pvector)))
    obspoint <- -log10(sort(pvector, decreasing=FALSE))

    eodata <- data.table(expect=expected, observe=obspoint)
    limax <- max(eodata)

    if(fullcut!=0){
      eodata <- eodata[expect > fullcut,]
    }

    size_point <- 0.1
    size_lamb <- 1.5
    modsize <- 5

    # Full plot
    plothole <- ggplot(eodata, aes(x = expect, y=observe))+
      xlab(expression(Expected~~-log[10](italic(p))))+
      ylab(expression(Observed~~-log[10](italic(p))))+
      annotate('segment', x = 0, y = 0, xend = limax, yend = limax, linewidth=0.2)+
      {if(cutoff !=0) geom_hline(yintercept = cutoff, color='red', linewidth=0.2)}+
      annotate('text', x=5, y=1.5, size= size_lamb, label=paste0("lambda %~% ~ ",smallGIF), parse=TRUE)+
      geom_point(size=size_point) +
      labs(title=fulltitle)+
      theme_minimal()+
      theme(plot.title = element_text(hjust = 0.5, size=modsize),
            axis.title = element_text(size=modsize - 1),
            axis.text = element_text(size=modsize / 1.5))

    return(plothole)

    # ggtitle(fulltitle, subtitle=bquote(lambda ~ " = " ~ .(GIF)))+

  }


  # setup
  mindful <- data.table(merged=c('DP_smoking_now', 'DP_poultry', 'DP_cheese', 'DP_tea', 'DP_BFP', 'SP_TVtime', 'SP_sleep_d', 'SP_veg_cook', 'SP_pork', 'SP_tea', 'SP_coffee', 'SP_sleep_dev', 'PP_TVtime', 'PP_sleep_dev', 'PP_salt', 'PP_tea'),
                        cutoff=c(5.5, 5.2, 4.14, 4.9, 4.7, 5.5, 3.94, 4.7, 4.23, 4.155, 4.7, 4.9, 4.33, 3.5, 4.5, 5.5))

  envlist <- c('Townsend', 'act0_d', 'TVtime', 'sleep_d','smoking_now', 'veg_cook', 'fish_oily', 'fish_lean', 'meat_proc', 'poultry', 'beef', 'lamb', 'pork', 'cheese', 'salt', 'tea', 'alc1', 'waist', 'getup', 'coffee', 'smoked_past', 'BFP', 'sleep_dev')

  allH <- 8
  allrow <- 6
  allcol <- 4
  sq <- allH / allrow
  allW <- sq * allcol

  allunits <- 'in'
  alldev <- 'png'

  # DP ----
  trait <- t <- 'DP'
  gg <- vector(mode='list', length=length(envlist))

  interject <- 0
  for(e in envlist){
    interject <- interject + 1
    print(paste0('Environment: ',e,' -------------------------------------'))
    gg[[interject]] <- claude(trait=t, env=e)
  }

  megadp <- ggarrange(plotlist=gg,
                      nrow=allrow,
                      ncol=allcol)
  # SP ----
  trait <- t <- 'SP'
  gg <- vector(mode='list', length=length(envlist))

  interject <- 0
  for(e in envlist){
    interject <- interject + 1
    print(paste0('Environment: ',e,' -------------------------------------'))
    gg[[interject]] <- claude(trait=t, env=e)
  }

  megasp <- ggarrange(plotlist=gg,
                      nrow=allrow,
                      ncol=allcol)
  # PP ----
  trait <- t <- 'PP'
  gg <- vector(mode='list', length=length(envlist))

  interject <- 0
  for(e in envlist){
    interject <- interject + 1
    print(paste0('Environment: ',e,' -------------------------------------'))
    gg[[interject]] <- claude(trait=t, env=e)
  }

  megapp <- ggarrange(plotlist=gg,
                      nrow=allrow,
                      ncol=allcol)


  # alldev   <- 'png'
  # outfile1 <- 'output/vault/sups/Figure_S3.png'
  # outfile2 <- 'output/vault/sups/Figure_S4.png'
  # outfile3 <- 'output/vault/sups/Figure_S5.png'
  outfile1 <- 'output/figs/final/Figure_S3.png'
  outfile2 <- 'output/figs/final/Figure_S4.png'
  outfile3 <- 'output/figs/final/Figure_S5.png'
  print('saving')
  ggsave(megadp, filename=outfile1, height=allH, width=allW, units=allunits, device=alldev, dpi=300)
  ggsave(megasp, filename=outfile2, height=allH, width=allW, units=allunits, device=alldev, dpi=300)
  ggsave(megapp, filename=outfile3, height=allH, width=allW, units=allunits, device=alldev, dpi=300)

}



print('Done!')
