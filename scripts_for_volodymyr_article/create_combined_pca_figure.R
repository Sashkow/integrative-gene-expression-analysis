#!/usr/bin/env Rscript

#' Create Combined PCA Figure with Proper Labels
#'
#' Regenerates 4 PCA plots with correct labels and combines into 2x2 figure
#' A - Before ComBat, Dataset
#' B - After ComBat, Dataset
#' C - Before ComBat, Trimester
#' D - After ComBat, Trimester
#'
#' @author Expression Integration Pipeline
#' @date 2025-11-15

cat("\n=== Creating Combined PCA Figure ===\n\n")

# Load required libraries
library(ggplot2)
library(gridExtra)
library(grid)

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

# Add panel labels A, B, C, D
add_panel_label <- function(plot, label) {
  plot +
    annotation_custom(
      grob = textGrob(label, x = 0.05, y = 0.95,
                      gp = gpar(fontsize = 16, fontface = "bold")),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
    )
}

plot_a <- add_panel_label(plot_a, "A")
plot_b <- add_panel_label(plot_b, "B")
plot_c <- add_panel_label(plot_c, "C")
plot_d <- add_panel_label(plot_d, "D")

# Combine into 2x2 grid
cat("\nCombining into 2x2 figure...\n")

# A4 half page dimensions: 8.27 x 5.85 inches
width_inches <- 8.27
height_inches <- 5.85

# Save as PNG
output_png <- file.path(output_dir, "combined_pca_figure.png")
cat("  Saving PNG:", output_png, "\n")
png(output_png, width = width_inches, height = height_inches, units = "in", res = 300)
grid.arrange(plot_a, plot_b, plot_c, plot_d, ncol = 2, nrow = 2)
dev.off()

# Save as PDF
output_pdf <- file.path(output_dir, "combined_pca_figure.pdf")
cat("  Saving PDF:", output_pdf, "\n")
pdf(output_pdf, width = width_inches, height = height_inches)
grid.arrange(plot_a, plot_b, plot_c, plot_d, ncol = 2, nrow = 2)
dev.off()

# Save as SVG
output_svg <- file.path(output_dir, "combined_pca_figure.svg")
cat("  Saving SVG:", output_svg, "\n")
svg(output_svg, width = width_inches, height = height_inches)
grid.arrange(plot_a, plot_b, plot_c, plot_d, ncol = 2, nrow = 2)
dev.off()

cat("\n✓ Combined figure created successfully\n")
cat("\nLayout:\n")
cat("  A (top-left):     Before ComBat - Dataset\n")
cat("  B (top-right):    After ComBat - Dataset\n")
cat("  C (bottom-left):  Before ComBat - Trimester\n")
cat("  D (bottom-right): After ComBat - Trimester\n\n")

cat("==================================================================\n")
cat("                    DONE                                          \n")
cat("==================================================================\n\n")
