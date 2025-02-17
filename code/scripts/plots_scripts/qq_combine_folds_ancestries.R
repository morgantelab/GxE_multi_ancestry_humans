rm(list=ls()); gc()
library(gridExtra)  # For arranging plots
library(ggpubr)     # For combining plots
library(ggplot2)    # For visualization
library(png)        # For reading PNG images
library(grid)       # For handling raster images

# Define directories
qqplot_dir <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/plots/"

# List of environments, BP traits, ancestries, and folds
envs <- c("Townsend", "act0_d", "TVtime", "sleep_d", "smoking_now", 
          "veg_cook", "fish_oily", "fish_lean", "meat_proc", "poultry", 
          "beef", "lamb", "pork", "cheese", "salt", "tea", "alc1", 
          "waist", "getup", "coffee", "smoked_past", "BFP", "sleep_dev")

bp_traits <- c("DP", "PP", "SP")  # BP traits only

ancestries <- c("mixed", "asian", "white", "black", "chinese")
folds <- 1:5

# Function to load PNG images as plots
load_qq_plot <- function(image_path) {
  if (file.exists(image_path)) {
    img <- readPNG(image_path)  # Read PNG
    return(rasterGrob(img, interpolate=TRUE))  # Convert to raster image for plotting
  } else {
    return(NULL)  # Return NULL if file is missing
  }
}

# ================== COMBINE FOLDS & ANCESTRIES FOR EACH BP TRAIT & ENV ================== #
for (trait in bp_traits) {
  for (env in envs) {
    print(paste("Grouping folds & ancestries for:", trait, "-", env))
    
    # Load fold QQ plots (row 1)
    fold_files <- sapply(folds, function(fold) {
      file.path(qqplot_dir, paste0("qqplot_gxe_", env, "_", fold, ".", trait, "0s.png"))
    })
    fold_plots <- lapply(fold_files, load_qq_plot)
    fold_plots <- fold_plots[!sapply(fold_plots, is.null)]  # Remove NULLs
    
    # Load ancestry QQ plots (row 2)
    ancestry_files <- sapply(ancestries, function(anc) {
      file.path(qqplot_dir, paste0("qqplot_gxe_", env, "_", anc, ".", trait, "0s.png"))
    })
    ancestry_plots <- lapply(ancestry_files, load_qq_plot)
    ancestry_plots <- ancestry_plots[!sapply(ancestry_plots, is.null)]  # Remove NULLs
    
    # Ensure both sets contain 5 plots before proceeding
    if (length(fold_plots) == 5 & length(ancestry_plots) == 5) {
      combined_plot <- ggarrange(plotlist=c(fold_plots, ancestry_plots), 
                                 ncol=5, nrow=2, 
                                 labels=c(paste("Fold", folds), ancestries))
      
      # Save the combined plot
      output_file <- file.path(qqplot_dir, paste0("qqplot_grouped_folds_ancestries_", trait, "_", env, ".png"))
      ggsave(output_file, combined_plot, width=15, height=8)
    }
  }
}

print("Grouped QQ plots (Folds & Ancestries) saved successfully!")
