#!/usr/bin/env Rscript

#' PCA Analysis: Before vs After ComBat Comparison
#'
#' Generates PCA plots for expression matrices before and after ComBat correction
#' For the GSE55439_GSE93520_GSE28551_GSE100051 dataset combination
#'
#' @author Expression Integration Pipeline
#' @date 2025-11-13

cat("\n=== PCA Analysis: Before vs After ComBat ===\n\n")

# Load required libraries
library(ggplot2)
library(factoextra)

# Source PCA analysis module
source("R/03_pca_analysis.R")

#' Custom PCA plot function using ggplot2 directly
#'
#' @param pca PCA object
#' @param pdata Phenodata
#' @param color_var Variable to color by
#' @param axes Principal components to plot
#' @param output_file Output file path (SVG)
#' @param output_dir_png PNG output directory
#' @param xlim Optional x-axis limits
#' @param ylim Optional y-axis limits
plot_pca_custom <- function(pca, pdata, color_var, axes = c(1, 2), output_file, output_dir_png = NULL, xlim = NULL, ylim = NULL) {

  # Extract PC scores
  pc_data <- as.data.frame(pca$x[, axes])
  colnames(pc_data) <- c("PC1", "PC2")

  # Add color variable
  pc_data$group <- as.factor(pdata[[color_var]])

  # Calculate variance explained
  var_explained <- summary(pca)$importance[2, axes] * 100

  # Create plot
  p <- ggplot(pc_data, aes(x = PC1, y = PC2, color = group)) +
    geom_point(size = 3, alpha = 0.7) +
    stat_ellipse(aes(group = group), type = "norm", level = 0.95, linetype = 2) +
    labs(
      title = paste("PCA colored by", color_var),
      x = sprintf("PC%d (%.1f%%)", axes[1], var_explained[1]),
      y = sprintf("PC%d (%.1f%%)", axes[2], var_explained[2]),
      color = color_var
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )

  # Apply fixed axis limits if provided
  if (!is.null(xlim)) {
    p <- p + xlim(xlim[1], xlim[2])
  }
  if (!is.null(ylim)) {
    p <- p + ylim(ylim[1], ylim[2])
  }

  # Save plot as SVG
  svg(output_file, width = 10, height = 8)
  print(p)
  dev.off()

  # Save plot as PNG in separate directory
  if (!is.null(output_dir_png)) {
    png_file <- file.path(output_dir_png, sub("\\.svg$", ".png", basename(output_file)))
  } else {
    png_file <- sub("\\.svg$", ".png", output_file)
  }
  png(png_file, width = 10, height = 8, units = "in", res = 300)
  print(p)
  dev.off()

  cat("Saved:", output_file, "and", png_file, "\n")
}

#' PCA plot with gestational age gradient
#'
#' @param pca PCA object
#' @param pdata Phenodata
#' @param axes Principal components to plot
#' @param output_file Output file path (SVG)
#' @param output_dir_png PNG output directory
#' @param xlim Optional x-axis limits
#' @param ylim Optional y-axis limits
plot_pca_ga_gradient <- function(pca, pdata, axes = c(1, 2), output_file, output_dir_png = NULL, xlim = NULL, ylim = NULL) {

  # Extract PC scores
  pc_data <- as.data.frame(pca$x[, axes])
  colnames(pc_data) <- c("PC1", "PC2")

  # Extract numeric gestational age
  ga_numeric <- as.numeric(gsub('[^0-9.]', '', pdata$Gestational.Age))
  pc_data$GA_weeks <- ga_numeric
  pc_data$GA_known <- !is.na(ga_numeric)
  pc_data$Trimester <- pdata$Gestational.Age.Category

  # Calculate variance explained
  var_explained <- summary(pca)$importance[2, axes] * 100

  # Define colors for gradient (same as trimester colors)
  # First trimester typically gets one color, second trimester another
  # We'll use a gradient from first trim color to second trim color
  color_first <- "#0073C2FF"  # Blue (typical for first group in jco palette)
  color_second <- "#EFC000FF"  # Orange/yellow (typical for second group)

  # Create plot
  p <- ggplot(pc_data, aes(x = PC1, y = PC2)) +
    # Plot known GA values with gradient
    geom_point(data = pc_data[pc_data$GA_known, ],
               aes(color = GA_weeks), size = 3, alpha = 0.8) +
    # Plot unknown GA values in black
    geom_point(data = pc_data[!pc_data$GA_known, ],
               color = "black", size = 3, alpha = 0.6, shape = 1) +
    scale_color_gradient(low = color_first, high = color_second,
                         name = "Gestational\nAge (weeks)",
                         na.value = "black") +
    labs(
      title = "PCA colored by Gestational Age (weeks)",
      subtitle = paste0("Known: ", sum(pc_data$GA_known), " samples, Unknown (black circles): ", sum(!pc_data$GA_known), " samples"),
      x = sprintf("PC%d (%.1f%%)", axes[1], var_explained[1]),
      y = sprintf("PC%d (%.1f%%)", axes[2], var_explained[2])
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )

  # Apply fixed axis limits if provided
  if (!is.null(xlim)) {
    p <- p + xlim(xlim[1], xlim[2])
  }
  if (!is.null(ylim)) {
    p <- p + ylim(ylim[1], ylim[2])
  }

  # Save plot as SVG
  svg(output_file, width = 10, height = 8)
  print(p)
  dev.off()

  # Save plot as PNG in separate directory
  if (!is.null(output_dir_png)) {
    png_file <- file.path(output_dir_png, sub("\\.svg$", ".png", basename(output_file)))
  } else {
    png_file <- sub("\\.svg$", ".png", output_file)
  }
  png(png_file, width = 10, height = 8, units = "in", res = 300)
  print(p)
  dev.off()

  cat("Saved:", output_file, "and", png_file, "\n")
}

# Define paths
data_dir <- "output/test_first_datasets/GSE55439_GSE93520_GSE28551_GSE100051"
pca_output_dir <- file.path(data_dir, "pca_plots")
pca_output_dir_png <- file.path(data_dir, "pca_plots_png")

# Create output directories
dir.create(pca_output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(pca_output_dir_png, showWarnings = FALSE, recursive = TRUE)
cat("Output directories:\n")
cat("  SVG:", pca_output_dir, "\n")
cat("  PNG:", pca_output_dir_png, "\n\n")

# Load phenodata
cat("Loading phenodata...\n")
pdata <- read.csv(file.path(data_dir, "phenodata.csv"), stringsAsFactors = FALSE)
cat("Samples:", nrow(pdata), "\n")

# Check available grouping variables
cat("\nGrouping variables:\n")
cat("  secondaryaccession:", length(unique(pdata$secondaryaccession)), "levels\n")
cat("  Gestational.Age.Category:", length(unique(pdata$Gestational.Age.Category)), "levels\n")
cat("  Combined.Fetus.Sex:", sum(!is.na(pdata$Combined.Fetus.Sex)), "non-NA values\n")

cat("\nDatasets in this analysis:\n")
print(table(pdata$secondaryaccession))

# Define variables to color by
color_vars <- c("secondaryaccession", "Gestational.Age.Category", "Combined.Fetus.Sex")

# ===================================================================
# ANALYSIS 1: BEFORE COMBAT
# ===================================================================

# Check if before/after combat files exist
before_combat_file <- file.path(data_dir, "merged_exprs_before_combat.tsv")
after_combat_file <- file.path(data_dir, "merged_exprs_after_combat.tsv")
old_format_file <- file.path(data_dir, "merged_exprs.tsv")

has_before_after <- file.exists(before_combat_file) && file.exists(after_combat_file)

if (!has_before_after) {
  if (file.exists(old_format_file)) {
    cat("\n⚠ Note: Using old file format (merged_exprs.tsv = after ComBat only)\n")
    cat("  To generate before/after files, re-run test_first_datasets.R\n\n")
    # Use old file as after-combat only
    after_combat_file <- old_format_file
    before_combat_file <- NULL
  } else {
    stop("No expression matrix files found in ", data_dir)
  }
}

if (!is.null(before_combat_file) && file.exists(before_combat_file)) {
  cat("\n")
  cat("==================================================================\n")
  cat("           PCA ANALYSIS: BEFORE COMBAT                           \n")
  cat("==================================================================\n\n")

  # Load expression matrix before ComBat
  cat("Loading expression matrix (before ComBat)...\n")
  exprs_before <- read.table(
    before_combat_file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE
  )
  cat("Expression matrix:", nrow(exprs_before), "genes x", ncol(exprs_before), "samples\n")

  # Align phenodata with expression matrix
  # Try different column names for sample matching
  sample_col <- if ("arraydatafile_exprscolumnnames" %in% colnames(pdata)) {
    "arraydatafile_exprscolumnnames"
  } else if ("exprs_column_names" %in% colnames(pdata)) {
    "exprs_column_names"
  } else {
    "Array.Data.File"
  }

  cat("Using sample column:", sample_col, "\n")
  pdata_aligned <- pdata[match(colnames(exprs_before), make.names(pdata[[sample_col]])), ]

  # Verify alignment
  if (any(is.na(pdata_aligned[[sample_col]]))) {
    stop("Sample alignment failed - some samples not found in phenodata")
  }
  cat("Phenodata aligned with expression matrix\n")

  # Perform PCA
  cat("\nPerforming PCA (before ComBat)...\n")
  pca_before <- perform_pca(exprs_before, center = TRUE, scale = FALSE)

  # Create scree plot
  cat("\nCreating scree plot...\n")
  plot_pca_scree(
    pca = pca_before,
    output_dir = pca_output_dir,
    output_dir_png = pca_output_dir_png,
    prefix = "before_combat",
    n_pcs = 20
  )

  cat("\n✓ Before ComBat analysis complete\n")
} else {
  cat("\n⚠ Skipping before-ComBat analysis (file not found)\n")
  pdata_aligned <- NULL
}

# ===================================================================
# ANALYSIS 2: AFTER COMBAT
# ===================================================================

cat("\n")
cat("==================================================================\n")
cat("           PCA ANALYSIS: AFTER COMBAT                            \n")
cat("==================================================================\n\n")

# Load expression matrix after ComBat
cat("Loading expression matrix (after ComBat)...\n")
exprs_after <- read.table(
  after_combat_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)
cat("Expression matrix:", nrow(exprs_after), "genes x", ncol(exprs_after), "samples\n")

# Align phenodata if not already done
if (is.null(pdata_aligned)) {
  # Try different column names for sample matching
  sample_col <- if ("arraydatafile_exprscolumnnames" %in% colnames(pdata)) {
    "arraydatafile_exprscolumnnames"
  } else if ("exprs_column_names" %in% colnames(pdata)) {
    "exprs_column_names"
  } else {
    "Array.Data.File"
  }

  cat("Using sample column:", sample_col, "\n")
  pdata_aligned <- pdata[match(colnames(exprs_after), make.names(pdata[[sample_col]])), ]
  if (any(is.na(pdata_aligned[[sample_col]]))) {
    stop("Sample alignment failed - some samples not found in phenodata")
  }
  cat("Phenodata aligned with expression matrix\n")
}

# Perform PCA
cat("\nPerforming PCA (after ComBat)...\n")
pca_after <- perform_pca(exprs_after, center = TRUE, scale = FALSE)

# Create scree plot
cat("\nCreating scree plot...\n")
plot_pca_scree(
  pca = pca_after,
  output_dir = pca_output_dir,
  output_dir_png = pca_output_dir_png,
  prefix = "after_combat",
  n_pcs = 20
)

cat("\n✓ After ComBat analysis complete\n")

# ===================================================================
# GENERATE PLOTS WITH OPTIMIZED AXIS SCALES
# ===================================================================

cat("\n")
cat("==================================================================\n")
cat("           GENERATING PLOTS WITH OPTIMIZED SCALES                \n")
cat("==================================================================\n\n")

if (exists("pca_before") && exists("pca_after")) {

  # Calculate separate axis limits for before-ComBat (PC1 vs PC2)
  cat("Calculating axis limits for before-ComBat (PC1 vs PC2)...\n")
  pc1_before <- range(pca_before$x[, 1])
  pc2_before <- range(pca_before$x[, 2])

  xlim_pc12_before <- pc1_before + c(-0.05, 0.05) * diff(pc1_before)
  ylim_pc12_before <- pc2_before + c(-0.05, 0.05) * diff(pc2_before)

  cat("  Before - PC1 limits:", round(xlim_pc12_before, 1), "\n")
  cat("  Before - PC2 limits:", round(ylim_pc12_before, 1), "\n")

  # Calculate separate axis limits for after-ComBat (PC1 vs PC2)
  cat("Calculating axis limits for after-ComBat (PC1 vs PC2)...\n")
  pc1_after <- range(pca_after$x[, 1])
  pc2_after <- range(pca_after$x[, 2])

  xlim_pc12_after <- pc1_after + c(-0.05, 0.05) * diff(pc1_after)
  ylim_pc12_after <- pc2_after + c(-0.05, 0.05) * diff(pc2_after)

  cat("  After - PC1 limits:", round(xlim_pc12_after, 1), "\n")
  cat("  After - PC2 limits:", round(ylim_pc12_after, 1), "\n")

  # Calculate separate axis limits for before-ComBat (PC3 vs PC4)
  cat("Calculating axis limits for before-ComBat (PC3 vs PC4)...\n")
  pc3_before <- range(pca_before$x[, 3])
  pc4_before <- range(pca_before$x[, 4])

  xlim_pc34_before <- pc3_before + c(-0.05, 0.05) * diff(pc3_before)
  ylim_pc34_before <- pc4_before + c(-0.05, 0.05) * diff(pc4_before)

  cat("  Before - PC3 limits:", round(xlim_pc34_before, 1), "\n")
  cat("  Before - PC4 limits:", round(ylim_pc34_before, 1), "\n")

  # Calculate separate axis limits for after-ComBat (PC3 vs PC4)
  cat("Calculating axis limits for after-ComBat (PC3 vs PC4)...\n")
  pc3_after <- range(pca_after$x[, 3])
  pc4_after <- range(pca_after$x[, 4])

  xlim_pc34_after <- pc3_after + c(-0.05, 0.05) * diff(pc3_after)
  ylim_pc34_after <- pc4_after + c(-0.05, 0.05) * diff(pc4_after)

  cat("  After - PC3 limits:", round(xlim_pc34_after, 1), "\n")
  cat("  After - PC4 limits:", round(ylim_pc34_after, 1), "\n\n")

  # Generate before-ComBat plots with appropriate scales
  cat("Creating before-ComBat plots (PC1 vs PC2)...\n")
  for (var in color_vars) {
    filename <- file.path(pca_output_dir, paste0("before_combat_", var, "_PC1_PC2.svg"))
    plot_pca_custom(pca_before, pdata_aligned, var, axes = c(1, 2),
                    output_file = filename, output_dir_png = pca_output_dir_png,
                    xlim = xlim_pc12_before, ylim = ylim_pc12_before)
  }

  cat("Creating before-ComBat plot (PC3 vs PC4)...\n")
  filename <- file.path(pca_output_dir, "before_combat_Combined.Fetus.Sex_PC3_PC4.svg")
  plot_pca_custom(pca_before, pdata_aligned, "Combined.Fetus.Sex", axes = c(3, 4),
                  output_file = filename, output_dir_png = pca_output_dir_png,
                  xlim = xlim_pc34_before, ylim = ylim_pc34_before)

  # Generate after-ComBat plots with appropriate scales
  cat("\nCreating after-ComBat plots (PC1 vs PC2)...\n")
  for (var in color_vars) {
    filename <- file.path(pca_output_dir, paste0("after_combat_", var, "_PC1_PC2.svg"))
    plot_pca_custom(pca_after, pdata_aligned, var, axes = c(1, 2),
                    output_file = filename, output_dir_png = pca_output_dir_png,
                    xlim = xlim_pc12_after, ylim = ylim_pc12_after)
  }

  cat("Creating after-ComBat plot (PC3 vs PC4)...\n")
  filename <- file.path(pca_output_dir, "after_combat_Combined.Fetus.Sex_PC3_PC4.svg")
  plot_pca_custom(pca_after, pdata_aligned, "Combined.Fetus.Sex", axes = c(3, 4),
                  output_file = filename, output_dir_png = pca_output_dir_png,
                  xlim = xlim_pc34_after, ylim = ylim_pc34_after)

  # Generate gestational age gradient plot
  cat("\nCreating after-ComBat plot with gestational age gradient (PC1 vs PC2)...\n")
  filename <- file.path(pca_output_dir, "after_combat_Gestational.Age_gradient_PC1_PC2.svg")
  plot_pca_ga_gradient(pca_after, pdata_aligned, axes = c(1, 2),
                       output_file = filename, output_dir_png = pca_output_dir_png,
                       xlim = xlim_pc12_after, ylim = ylim_pc12_after)

  cat("\n✓ All plots generated with appropriate scales\n")

} else if (exists("pca_after")) {
  # Only after-ComBat available, generate plots without shared scales
  cat("Generating after-ComBat plots (no before-ComBat data for comparison)...\n")
  for (var in color_vars) {
    filename <- file.path(pca_output_dir, paste0("after_combat_", var, "_PC1_PC2.svg"))
    plot_pca_custom(pca_after, pdata_aligned, var, axes = c(1, 2), output_file = filename)
  }
  filename <- file.path(pca_output_dir, "after_combat_Combined.Fetus.Sex_PC3_PC4.svg")
  plot_pca_custom(pca_after, pdata_aligned, "Combined.Fetus.Sex", axes = c(3, 4), output_file = filename)
}

# ===================================================================
# SUMMARY
# ===================================================================

cat("\n")
cat("==================================================================\n")
cat("                        SUMMARY                                   \n")
cat("==================================================================\n\n")

if (exists("pca_before")) {
  cat("Before ComBat - Variance explained by top 5 PCs:\n")
  var_before <- summary(pca_before)$importance[2, 1:5] * 100
  for (i in 1:5) {
    cat(sprintf("  PC%d: %.2f%%\n", i, var_before[i]))
  }
  cat("\n")
}

cat("After ComBat - Variance explained by top 5 PCs:\n")
var_after <- summary(pca_after)$importance[2, 1:5] * 100
for (i in 1:5) {
  cat(sprintf("  PC%d: %.2f%%\n", i, var_after[i]))
}

cat("\nAll plots saved to:", pca_output_dir, "\n")

cat("\nPCA plots generated:\n")
if (exists("pca_before")) {
  cat("  Before ComBat:\n")
  cat("    - before_combat_scree.svg (scree plot)\n")
  cat("    - before_combat_secondaryaccession_PC1_PC2.svg (colored by dataset)\n")
  cat("    - before_combat_Gestational.Age.Category_PC1_PC2.svg (colored by trimester)\n")
  cat("    - before_combat_Combined.Fetus.Sex_PC1_PC2.svg (colored by fetus sex, PC1 vs PC2)\n")
  cat("    - before_combat_Combined.Fetus.Sex_PC3_PC4.svg (colored by fetus sex, PC3 vs PC4)\n")
  cat("\n")
}
cat("  After ComBat:\n")
cat("    - after_combat_scree.svg (scree plot)\n")
cat("    - after_combat_secondaryaccession_PC1_PC2.svg (colored by dataset)\n")
cat("    - after_combat_Gestational.Age.Category_PC1_PC2.svg (colored by trimester)\n")
cat("    - after_combat_Gestational.Age_gradient_PC1_PC2.svg (colored by exact gestational age with gradient)\n")
cat("    - after_combat_Combined.Fetus.Sex_PC1_PC2.svg (colored by fetus sex, PC1 vs PC2)\n")
cat("    - after_combat_Combined.Fetus.Sex_PC3_PC4.svg (colored by fetus sex, PC3 vs PC4)\n")

cat("\n")
cat("==================================================================\n")
cat("                    ANALYSIS COMPLETE                             \n")
cat("==================================================================\n\n")
