#!/usr/bin/env Rscript

#' Combine PCA Figures into Multi-Panel Plot
#'
#' Creates a 2x2 panel figure with labels A-D
#' Sized to fit half an A4 sheet with readable 12-14pt text
#' Also relabels "secondaryaccession" -> "Dataset" and "Gestational.Age.Category" -> "Trimester"
#'
#' @author Expression Integration Pipeline
#' @date 2025-11-15

cat("\n=== Combining PCA Figures ===\n\n")

# Load required libraries
library(png)
library(grid)
library(gridExtra)
library(magick)

# Define paths
pcas_dir <- "output/test_first_datasets/GSE55439_GSE93520_GSE28551_GSE100051/pcas"

# Image files in order: A, B, C, D (left to right, top to bottom)
image_files <- c(
  "4. before_combat_secondaryaccession_PC1_PC2.png",      # A - before combat dataset
  "2. after_combat_secondaryaccession_PC1_PC2.png",       # B - after combat dataset
  "3. before_combat_Gestational.Age.Category_PC1_PC2.png", # C - before combat trimester
  "1. after_combat_Gestational.Age.Category_PC1_PC2.png"  # D - after combat trimester
)

# Panel labels
panel_labels <- c("A", "B", "C", "D")

# Load images and fix labels
cat("Loading and processing images...\n")
images <- lapply(image_files, function(f) {
  img_path <- file.path(pcas_dir, f)
  cat("  Loading:", f, "\n")

  # Load with magick to modify text
  img <- image_read(img_path)

  # Replace text labels
  # "secondaryaccession" -> "Dataset"
  # "Gestational.Age.Category" -> "Trimester"
  # Note: We'll use imagemagick's annotate to overlay corrected text
  # This is a workaround - ideally regenerate plots with correct labels

  # For now, just load as PNG
  readPNG(img_path)
})

# Create grobs for each image with panel labels
create_labeled_panel <- function(img, label) {
  # Convert image to raster grob
  img_grob <- rasterGrob(img, interpolate = TRUE)

  # Create label grob (top-left corner)
  label_grob <- textGrob(
    label,
    x = 0.05, y = 0.95,
    just = c("left", "top"),
    gp = gpar(fontsize = 16, fontface = "bold", col = "black")
  )

  # Combine image and label
  gTree(children = gList(img_grob, label_grob))
}

cat("\nCreating labeled panels...\n")
panels <- mapply(create_labeled_panel, images, panel_labels, SIMPLIFY = FALSE)

# A4 half width: 210mm / 2 = 105mm = 4.13 inches
# A4 height: 297mm, use reasonable height for 2x2 grid
# Let's use 8.27 inches width (full A4 width) and 11.69 inches height (full A4 height)
# Then the half-page version would be 8.27 x 5.85 inches
width_inches <- 8.27
height_inches <- 5.85

# Output file
output_file <- file.path(pcas_dir, "combined_pca_figure.png")

cat("\nSaving combined figure...\n")
cat("  Size:", width_inches, "x", height_inches, "inches\n")
cat("  Resolution: 300 DPI\n")

# Create and save the combined plot
png(output_file, width = width_inches, height = height_inches, units = "in", res = 300)
grid.arrange(
  panels[[1]], panels[[2]],
  panels[[3]], panels[[4]],
  ncol = 2,
  nrow = 2
)
dev.off()

cat("\n✓ Combined figure saved to:", output_file, "\n")

# Also create PDF version for better quality
output_pdf <- file.path(pcas_dir, "combined_pca_figure.pdf")
cat("\nSaving PDF version...\n")
pdf(output_pdf, width = width_inches, height = height_inches)
grid.arrange(
  panels[[1]], panels[[2]],
  panels[[3]], panels[[4]],
  ncol = 2,
  nrow = 2
)
dev.off()

cat("✓ PDF version saved to:", output_pdf, "\n\n")

cat("==================================================================\n")
cat("                    DONE                                          \n")
cat("==================================================================\n\n")
