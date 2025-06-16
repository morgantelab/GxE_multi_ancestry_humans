rm(list=ls()); gc()
library(gridExtra)  # For combining plots
library(ggpubr)     # For arranging plots
library(ggplot2)    # For visualization
library(png)        # For reading PNG images
library(grid)       # For handling raster images

# Define directories
qqplot_dir <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/plots/"

# List of environments, traits, and ancestries
envs <- c("Townsend", "act0_d", "TVtime", "sleep_d", "smoking_now", 
          "veg_cook", "fish_oily", "fish_lean", "meat_proc", "poultry", 
          "beef", "lamb", "pork", "cheese", "salt", "tea", "alc1", 
          "waist", "getup", "coffee", "smoked_past", "BFP", "sleep_dev")

traits <- c("DP", "PP", "SP")

ancestries <- c("mixed", "asian", "white", "black", "chinese")

# Function to load and display PNG images
load_qq_plot <- function(image_path) {
  img <- readPNG(image_path)  # Read PNG
  rasterGrob(img, interpolate=TRUE)  # Convert to a raster image for plotting
}

# ============ Group 5 Ancestries for the Same Trait & Environment ============ #
for (trait in traits) {
  for (env in envs) {
    print(paste("Grouping ancestries for:", trait, "-", env))
    
    # Find QQ plots for the 5 ancestries of the given trait & environment
    ancestry_files <- list.files(qqplot_dir, pattern=paste0("gxe_", env, "_(", paste(ancestries, collapse="|"), ")\\.", trait, "0s\\.png$"), full.names=TRUE)
    
    # Ensure we have all 5 ancestries
    ancestry_files <- ancestry_files[order(match(gsub(".*_(.*?)\\..*\\.png$", "\\1", ancestry_files), ancestries))]
    
    if (length(ancestry_files) == 5) {
      plot_list <- lapply(ancestry_files, load_qq_plot)  # Load images
      grouped_plot <- ggarrange(plotlist=plot_list, ncol=5, nrow=1, labels=ancestries)  # Arrange ancestries side-by-side
      
      # Save grouped image in the same directory
      output_file <- file.path(qqplot_dir, paste0("qqplot_grouped_ancestries_2", trait, "_", env, ".png"))
      ggsave(output_file, grouped_plot, width=15, height=5)
    } else {
      print(paste("Skipping:", trait, "-", env, "- Missing ancestries"))
    }
  }
}

print("Grouped ancestry QQ plots saved successfully!")
