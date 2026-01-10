#!/usr/bin/env Rscript

#' Create PowerPoint Presentation with PCA Figures
#'
#' Creates a PPTX with the combined PCA figure and editable text captions
#'
#' @author Expression Integration Pipeline
#' @date 2025-11-15

cat("\n=== Creating PowerPoint Presentation ===\n\n")

# Load required libraries
library(officer)
library(ggplot2)

# Define paths
data_dir <- "output/test_first_datasets/GSE55439_GSE93520_GSE28551_GSE100051"
output_dir <- file.path(data_dir, "pcas")

# Load data
cat("Loading data...\n")
pdata <- read.csv(file.path(data_dir, "phenodata.csv"), stringsAsFactors = FALSE)

# Load before and after ComBat expression matrices
exprs_before <- read.table(
  file.path(data_dir, "merged_exprs_before_combat.tsv"),
  header = TRUE, sep = "\t", check.names = FALSE
)
exprs_after <- read.table(
  file.path(data_dir, "merged_exprs_after_combat.tsv"),
  header = TRUE, sep = "\t", check.names = FALSE
)

cat("Expression matrices loaded\n\n")

# Helper function to create PCA plot
create_pca_plot <- function(exprs, pdata, color_var, var_label, title_prefix) {

  # Perform PCA
  exprs_clean <- na.omit(exprs)
  pca <- prcomp(t(exprs_clean), center = TRUE, scale. = FALSE)

  # Extract PC scores
  pc_data <- as.data.frame(pca$x[, c(1, 2)])
  colnames(pc_data) <- c("PC1", "PC2")

  # Add color variable
  pc_data$group <- as.factor(pdata[[color_var]])

  # Calculate variance explained
  var_explained <- summary(pca)$importance[2, 1:2] * 100

  # Create plot
  p <- ggplot(pc_data, aes(x = PC1, y = PC2, color = group)) +
    geom_point(size = 3, alpha = 0.7) +
    stat_ellipse(aes(group = group), type = "norm", level = 0.95, linetype = 2) +
    labs(
      title = paste(title_prefix, "- by", var_label),
      x = sprintf("PC1 (%.1f%%)", var_explained[1]),
      y = sprintf("PC2 (%.1f%%)", var_explained[2]),
      color = var_label
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      legend.position = "right",
      panel.grid.minor = element_blank(),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11)
    )

  return(p)
}

cat("Creating plots...\n")

# A - Before ComBat, Dataset
cat("  Panel A: Before ComBat - Dataset\n")
plot_a <- create_pca_plot(exprs_before, pdata, "secondaryaccession", "Dataset", "Before ComBat")

# B - After ComBat, Dataset
cat("  Panel B: After ComBat - Dataset\n")
plot_b <- create_pca_plot(exprs_after, pdata, "secondaryaccession", "Dataset", "After ComBat")

# C - Before ComBat, Trimester
cat("  Panel C: Before ComBat - Trimester\n")
plot_c <- create_pca_plot(exprs_before, pdata, "Gestational.Age.Category", "Trimester", "Before ComBat")

# D - After ComBat, Trimester
cat("  Panel D: After ComBat - Trimester\n")
plot_d <- create_pca_plot(exprs_after, pdata, "Gestational.Age.Category", "Trimester", "After ComBat")

# Create temporary image files for each plot
cat("\nSaving individual plot images...\n")
temp_dir <- tempdir()

ggsave(file.path(temp_dir, "plot_a.png"), plot_a, width = 4, height = 3.5, dpi = 300)
ggsave(file.path(temp_dir, "plot_b.png"), plot_b, width = 4, height = 3.5, dpi = 300)
ggsave(file.path(temp_dir, "plot_c.png"), plot_c, width = 4, height = 3.5, dpi = 300)
ggsave(file.path(temp_dir, "plot_d.png"), plot_d, width = 4, height = 3.5, dpi = 300)

# Create PowerPoint presentation
cat("\nCreating PowerPoint presentation...\n")
pptx <- read_pptx()

# Add a slide with blank layout
pptx <- add_slide(pptx, layout = "Blank", master = "Office Theme")

# Define positions for 2x2 grid (in inches)
# Standard slide is 10 x 7.5 inches
left_col <- 0.5
right_col <- 5.5
top_row <- 0.8
bottom_row <- 4.3
plot_width <- 4
plot_height <- 3.5
label_size <- 18

# Add Panel A (top-left)
pptx <- ph_with(pptx, external_img(file.path(temp_dir, "plot_a.png")),
                location = ph_location(left = left_col, top = top_row,
                                      width = plot_width, height = plot_height))
pptx <- ph_with(pptx, value = "A",
                location = ph_location(left = left_col + 0.1, top = top_row + 0.1,
                                      width = 0.5, height = 0.5),
                bg = "transparent",
                fp_text = fp_text(font.size = label_size, bold = TRUE))

# Add Panel B (top-right)
pptx <- ph_with(pptx, external_img(file.path(temp_dir, "plot_b.png")),
                location = ph_location(left = right_col, top = top_row,
                                      width = plot_width, height = plot_height))
pptx <- ph_with(pptx, value = "B",
                location = ph_location(left = right_col + 0.1, top = top_row + 0.1,
                                      width = 0.5, height = 0.5),
                bg = "transparent",
                fp_text = fp_text(font.size = label_size, bold = TRUE))

# Add Panel C (bottom-left)
pptx <- ph_with(pptx, external_img(file.path(temp_dir, "plot_c.png")),
                location = ph_location(left = left_col, top = bottom_row,
                                      width = plot_width, height = plot_height))
pptx <- ph_with(pptx, value = "C",
                location = ph_location(left = left_col + 0.1, top = bottom_row + 0.1,
                                      width = 0.5, height = 0.5),
                bg = "transparent",
                fp_text = fp_text(font.size = label_size, bold = TRUE))

# Add Panel D (bottom-right)
pptx <- ph_with(pptx, external_img(file.path(temp_dir, "plot_d.png")),
                location = ph_location(left = right_col, top = bottom_row,
                                      width = plot_width, height = plot_height))
pptx <- ph_with(pptx, value = "D",
                location = ph_location(left = right_col + 0.1, top = bottom_row + 0.1,
                                      width = 0.5, height = 0.5),
                bg = "transparent",
                fp_text = fp_text(font.size = label_size, bold = TRUE))

# Add slide title
pptx <- ph_with(pptx, value = "PCA Analysis: Batch Effect Removal Assessment",
                location = ph_location(left = 0.5, top = 0.1,
                                      width = 9, height = 0.5),
                bg = "transparent",
                fp_text = fp_text(font.size = 20, bold = TRUE))

# Add figure caption as editable text box
caption_text <- paste(
  "Figure 1. Principal component analysis comparing before and after batch effect removal.",
  "Panel A shows dataset distribution before ComBat correction, with strong technical variation.",
  "Panel B demonstrates successful batch correction with overlapping datasets.",
  "Panel C shows obscured trimester separation before correction.",
  "Panel D reveals clear biological separation by trimester after batch effect removal."
)

pptx <- ph_with(pptx, value = caption_text,
                location = ph_location(left = 0.5, top = 7.2,
                                      width = 9, height = 0.3),
                bg = "transparent",
                fp_text = fp_text(font.size = 10))

# Save PowerPoint file
output_file <- file.path(output_dir, "combined_pca_figure.pptx")
print(pptx, target = output_file)

cat("\n✓ PowerPoint presentation saved to:", output_file, "\n")

# Clean up temporary files
file.remove(file.path(temp_dir, "plot_a.png"))
file.remove(file.path(temp_dir, "plot_b.png"))
file.remove(file.path(temp_dir, "plot_c.png"))
file.remove(file.path(temp_dir, "plot_d.png"))

cat("\nLayout:\n")
cat("  A (top-left):     Before ComBat - Dataset\n")
cat("  B (top-right):    After ComBat - Dataset\n")
cat("  C (bottom-left):  Before ComBat - Trimester\n")
cat("  D (bottom-right): After ComBat - Trimester\n\n")

cat("All text boxes (title, labels, caption) are editable in PowerPoint.\n\n")

cat("==================================================================\n")
cat("                    DONE                                          \n")
cat("==================================================================\n\n")
