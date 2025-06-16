rm(list=ls()); gc()
library(gridExtra)  # For combining plots
library(ggpubr)     # For arranging plots
library(ggplot2)    # For visualization
library(png)        # For reading PNG images
library(grid)       # For handling raster images

# Define directories
qqplot_dir <- "/data2/morgante_lab/ukbiobank_projects/GxE_multi_ancestry/plots/"

# List of environments, BP traits, and ancestries
envs <- c("Townsend", "act0_d", "TVtime", "sleep_d", "smoking_now", 
          "veg_cook", "fish_oily", "fish_lean", "meat_proc", "poultry", 
          "beef", "lamb", "pork", "cheese", "salt", "tea", "alc1", 
          "waist", "getup", "coffee", "smoked_past", "BFP", "sleep_dev")

bp_traits <- c("DP", "PP", "SP")  # Only BP traits

ancestries <- c("mixed", "asian", "white", "black", "chinese")

# Function to load PNG images as plots
load_qq_plot <- function(image_path) {
  img <- readPNG(image_path)  # Read PNG
  rasterGrob(img, interpolate=TRUE)  # Convert to a raster image for plotting
}

# ============ STEP: Group 5 Ancestries for Each BP Trait & Environment ============ #
for (trait in bp_traits) {
  for (env in envs) {
    print(paste("Grouping ancestries for:", trait, "-", env))
    
    # Find all ancestry-specific QQ plots for this trait & environment
    ancestry_files <- sapply(ancestries, function(anc) {
      file.path(qqplot_dir, paste0("qqplot_gxe_", env, "_", anc, ".", trait, "0s.png"))
    })
    
    # Keep only existing files
    ancestry_files <- ancestry_files[file.exists(ancestry_files)]
    
    if (length(ancestry_files) == 5) {  # Ensure all 5 ancestries are present
      plot_list <- lapply(ancestry_files, load_qq_plot)  # Load images
      
      # Arrange in a 2x3 grid (or 1x5 if you prefer)
      grouped_plot <- ggarrange(plotlist=plot_list, ncol=3, nrow=2, labels=ancestries)
      
      # Save grouped image in the same directory
      output_file <- file.path(qqplot_dir, paste0("qqplot_grouped_ancestries_", trait, "_", env, ".png"))
      ggsave(output_file, grouped_plot, width=12, height=8)
    }
  }
}

print("Grouped QQ plots saved successfully!")
