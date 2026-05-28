#!/usr/bin/env Rscript

#' Run pipeline for 2nd vs 3rd Trimester (6 datasets) and create PCA multiplot
#'
#' Datasets: GSE37901, GSE6573, GSE9984, GSE73374, GSE73685 (Affymetrix), GSE35574 (Illumina)

cat("\n=== Step 1: Running Pipeline ===\n\n")

source("scripts/analysis/compare_baseline_addons.R")

run_dataset_comparison(
  config_file = "config/volodymyr/config_volodymyr_2_3_trim.yaml",
  output_dir = "output/volodymyr/volodymyr_2_3_trim",
  trimester_col_1 = "Second Trimester",
  trimester_col_2 = "Term",
  addon_config_key = "term_datasets",
  file_prefix = "term_",
  type_label = "with_term",
  enable_logging = FALSE,
  use_global_exclusions = FALSE
)

cat("\n=== Step 2: Creating Combined PCA Figure ===\n\n")

library(ggplot2)
library(gridExtra)
library(grid)

data_dir <- "output/volodymyr/volodymyr_2_3_trim/GSE35574"
output_dir <- file.path(data_dir, "pcas")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load data
cat("Loading data...\n")
pdata <- read.csv(file.path(data_dir, "phenodata.csv"), stringsAsFactors = FALSE)

exprs_before <- read.table(
  file.path(data_dir, "merged_exprs_before_combat.tsv"),
  header = TRUE, sep = "\t", check.names = FALSE
)
exprs_after <- read.table(
  file.path(data_dir, "merged_exprs_after_combat.tsv"),
  header = TRUE, sep = "\t", check.names = FALSE
)

cat("Samples:", nrow(pdata), "\n")
cat("Datasets:", paste(unique(pdata$secondaryaccession), collapse = ", "), "\n")
cat("Genes before:", nrow(exprs_before), " after:", nrow(exprs_after), "\n\n")

# PCA plot helper
create_pca_plot <- function(exprs, pdata, color_var, var_label, title_prefix) {
  exprs_clean <- na.omit(exprs)
  pca <- prcomp(t(exprs_clean), center = TRUE, scale. = FALSE)

  pc_data <- as.data.frame(pca$x[, c(1, 2)])
  colnames(pc_data) <- c("PC1", "PC2")
  pc_data$group <- as.factor(pdata[[color_var]])

  var_explained <- summary(pca)$importance[2, 1:2] * 100

  ggplot(pc_data, aes(x = PC1, y = PC2, color = group)) +
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
}

cat("Creating plots...\n")

plot_a <- create_pca_plot(exprs_before, pdata, "secondaryaccession", "Dataset", "Before ComBat")
plot_b <- create_pca_plot(exprs_after, pdata, "secondaryaccession", "Dataset", "After ComBat")
plot_c <- create_pca_plot(exprs_before, pdata, "Gestational.Age.Category", "Trimester", "Before ComBat")
plot_d <- create_pca_plot(exprs_after, pdata, "Gestational.Age.Category", "Trimester", "After ComBat")

# Add panel labels
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
width_inches <- 8.27
height_inches <- 5.85

output_png <- file.path(output_dir, "combined_pca_figure.png")
cat("  Saving PNG:", output_png, "\n")
png(output_png, width = width_inches, height = height_inches, units = "in", res = 300)
grid.arrange(plot_a, plot_b, plot_c, plot_d, ncol = 2, nrow = 2)
dev.off()

output_pdf <- file.path(output_dir, "combined_pca_figure.pdf")
cat("  Saving PDF:", output_pdf, "\n")
pdf(output_pdf, width = width_inches, height = height_inches)
grid.arrange(plot_a, plot_b, plot_c, plot_d, ncol = 2, nrow = 2)
dev.off()

output_svg <- file.path(output_dir, "combined_pca_figure.svg")
cat("  Saving SVG:", output_svg, "\n")
svg(output_svg, width = width_inches, height = height_inches)
grid.arrange(plot_a, plot_b, plot_c, plot_d, ncol = 2, nrow = 2)
dev.off()

cat("\n=== DONE ===\n")
cat("Layout:\n")
cat("  A (top-left):     Before ComBat - Dataset\n")
cat("  B (top-right):    After ComBat - Dataset\n")
cat("  C (bottom-left):  Before ComBat - Trimester\n")
cat("  D (bottom-right): After ComBat - Trimester\n\n")
cat("Output:", output_dir, "\n")
