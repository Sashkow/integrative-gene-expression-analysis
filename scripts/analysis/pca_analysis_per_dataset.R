#!/usr/bin/env Rscript

#' PCA Analysis: Individual Dataset Analysis
#'
#' Generates PCA plots for each dataset individually
#' Plots include: by trimester, by sex, by gestational age gradient, and by sample name
#'
#' @author Expression Integration Pipeline
#' @date 2025-11-13

cat("\n=== PCA Analysis: Individual Datasets ===\n\n")

# Load required libraries
library(ggplot2)
library(ggrepel)  # For non-overlapping labels

# Source PCA analysis module
source("R/03_pca_analysis.R")

#' PCA plot with sample labels
#'
#' @param pca PCA object
#' @param pdata Phenodata
#' @param sample_col Column containing sample names
#' @param axes Principal components to plot
#' @param output_file Output file path (SVG)
#' @param output_dir_png PNG output directory
#' @param xlim Optional x-axis limits
#' @param ylim Optional y-axis limits
plot_pca_samples <- function(pca, pdata, sample_col, axes = c(1, 2), output_file, output_dir_png = NULL, xlim = NULL, ylim = NULL) {

  # Extract PC scores
  pc_data <- as.data.frame(pca$x[, axes])
  colnames(pc_data) <- c("PC1", "PC2")

  # Add sample names
  pc_data$sample <- pdata[[sample_col]]

  # Calculate variance explained
  var_explained <- summary(pca)$importance[2, axes] * 100

  # Create plot with non-overlapping labels
  p <- ggplot(pc_data, aes(x = PC1, y = PC2, label = sample)) +
    geom_point(size = 2, alpha = 0.6, color = "steelblue") +
    geom_text_repel(
      size = 2.5,
      max.overlaps = Inf,
      box.padding = 0.5,
      point.padding = 0.3,
      segment.color = "grey50",
      segment.size = 0.2
    ) +
    labs(
      title = "PCA with Sample Labels",
      x = sprintf("PC%d (%.1f%%)", axes[1], var_explained[1]),
      y = sprintf("PC%d (%.1f%%)", axes[2], var_explained[2])
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
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
  svg(output_file, width = 12, height = 10)
  print(p)
  dev.off()

  # Save plot as PNG in separate directory
  if (!is.null(output_dir_png)) {
    png_file <- file.path(output_dir_png, sub("\\.svg$", ".png", basename(output_file)))
  } else {
    png_file <- sub("\\.svg$", ".png", output_file)
  }
  png(png_file, width = 12, height = 10, units = "in", res = 300)
  print(p)
  dev.off()

  cat("Saved:", output_file, "and", png_file, "\n")
}

#' Custom PCA plot function for individual datasets
plot_pca_custom_single <- function(pca, pdata, color_var, axes = c(1, 2), output_file, output_dir_png = NULL, xlim = NULL, ylim = NULL) {

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

#' PCA plot with gestational age gradient for individual datasets
plot_pca_ga_gradient_single <- function(pca, pdata, axes = c(1, 2), output_file, output_dir_png = NULL, xlim = NULL, ylim = NULL) {

  # Extract PC scores
  pc_data <- as.data.frame(pca$x[, axes])
  colnames(pc_data) <- c("PC1", "PC2")

  # Extract numeric gestational age
  ga_numeric <- as.numeric(gsub('[^0-9.]', '', pdata$Gestational.Age))
  pc_data$GA_weeks <- ga_numeric
  pc_data$GA_known <- !is.na(ga_numeric)

  # Calculate variance explained
  var_explained <- summary(pca)$importance[2, axes] * 100

  # Define colors for gradient
  color_first <- "#0073C2FF"  # Blue
  color_second <- "#EFC000FF"  # Orange/yellow

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
pca_output_base <- file.path(data_dir, "pca_plots_per_dataset")
pca_output_base_png <- file.path(data_dir, "pca_plots_per_dataset_png")

# Create base output directories
dir.create(pca_output_base, showWarnings = FALSE, recursive = TRUE)
dir.create(pca_output_base_png, showWarnings = FALSE, recursive = TRUE)
cat("Output directories:\n")
cat("  SVG:", pca_output_base, "\n")
cat("  PNG:", pca_output_base_png, "\n\n")

# Load phenodata and expression matrix
cat("Loading data...\n")
pdata <- read.csv(file.path(data_dir, "phenodata.csv"), stringsAsFactors = FALSE)
exprs <- read.table(
  file.path(data_dir, "merged_exprs_after_combat.tsv"),
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)
cat("Expression matrix:", nrow(exprs), "genes x", ncol(exprs), "samples\n\n")

# Get unique datasets
datasets <- unique(pdata$secondaryaccession)
datasets <- datasets[datasets != ""]  # Remove empty entries
cat("Found", length(datasets), "datasets:\n")
print(datasets)
cat("\n")

# Analyze each dataset
for (dataset_id in datasets) {

  cat("\n")
  cat("==================================================================\n")
  cat("           DATASET:", dataset_id, "\n")
  cat("==================================================================\n\n")

  # Filter phenodata for this dataset
  pdata_dataset <- pdata[pdata$secondaryaccession == dataset_id, ]
  n_samples <- nrow(pdata_dataset)

  cat("Samples:", n_samples, "\n")

  # Check if dataset has both trimesters
  trimesters <- unique(pdata_dataset$Gestational.Age.Category)
  n_trimesters <- length(trimesters[!is.na(trimesters) & trimesters != ""])

  cat("Trimesters:", paste(trimesters[!is.na(trimesters) & trimesters != ""], collapse = ", "), "\n")

  # Skip if too few samples
  if (n_samples < 3) {
    cat("⚠ Too few samples (", n_samples, ") - skipping\n")
    next
  }

  # Determine sample column
  sample_col <- if ("arraydatafile_exprscolumnnames" %in% colnames(pdata_dataset)) {
    "arraydatafile_exprscolumnnames"
  } else if ("exprs_column_names" %in% colnames(pdata_dataset)) {
    "exprs_column_names"
  } else {
    "Array.Data.File"
  }

  # Get expression data for this dataset
  sample_names <- make.names(pdata_dataset[[sample_col]])
  exprs_dataset <- exprs[, colnames(exprs) %in% sample_names]

  cat("Expression matrix for dataset:", nrow(exprs_dataset), "genes x", ncol(exprs_dataset), "samples\n")

  # Align phenodata
  pdata_aligned <- pdata_dataset[match(colnames(exprs_dataset), sample_names), ]

  # Create dataset-specific output directories
  dataset_dir <- file.path(pca_output_base, dataset_id)
  dataset_dir_png <- file.path(pca_output_base_png, dataset_id)
  dir.create(dataset_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(dataset_dir_png, showWarnings = FALSE, recursive = TRUE)

  # Perform PCA
  cat("\nPerforming PCA...\n")
  pca_dataset <- perform_pca(exprs_dataset, center = TRUE, scale = FALSE)

  # Calculate axis limits
  pc1_range <- range(pca_dataset$x[, 1])
  pc2_range <- range(pca_dataset$x[, 2])
  xlim <- pc1_range + c(-0.05, 0.05) * diff(pc1_range)
  ylim <- pc2_range + c(-0.05, 0.05) * diff(pc2_range)

  # Create scree plot
  cat("Creating scree plot...\n")
  plot_pca_scree(
    pca = pca_dataset,
    output_dir = dataset_dir,
    output_dir_png = dataset_dir_png,
    prefix = dataset_id,
    n_pcs = min(10, ncol(exprs_dataset) - 1)
  )

  # Plot by trimester (if multiple trimesters)
  if (n_trimesters >= 2) {
    cat("Creating PCA plot by trimester...\n")
    filename <- file.path(dataset_dir, paste0(dataset_id, "_Gestational.Age.Category_PC1_PC2.svg"))
    tryCatch({
      plot_pca_custom_single(pca_dataset, pdata_aligned, "Gestational.Age.Category",
                             axes = c(1, 2), output_file = filename,
                             output_dir_png = dataset_dir_png,
                             xlim = xlim, ylim = ylim)
    }, error = function(e) {
      cat("⚠ Error creating trimester plot:", conditionMessage(e), "\n")
    })
  } else {
    cat("⚠ Only one trimester - skipping trimester plot\n")
  }

  # Plot by sex
  n_sex <- sum(!is.na(pdata_aligned$Combined.Fetus.Sex) & pdata_aligned$Combined.Fetus.Sex != "")
  if (n_sex > 0) {
    cat("Creating PCA plot by sex...\n")
    filename <- file.path(dataset_dir, paste0(dataset_id, "_Combined.Fetus.Sex_PC1_PC2.svg"))
    tryCatch({
      plot_pca_custom_single(pca_dataset, pdata_aligned, "Combined.Fetus.Sex",
                             axes = c(1, 2), output_file = filename,
                             output_dir_png = dataset_dir_png,
                             xlim = xlim, ylim = ylim)
    }, error = function(e) {
      cat("⚠ Error creating sex plot:", conditionMessage(e), "\n")
    })
  } else {
    cat("⚠ No sex information - skipping sex plot\n")
  }

  # Plot by gestational age gradient
  ga_numeric <- as.numeric(gsub('[^0-9.]', '', pdata_aligned$Gestational.Age))
  n_ga <- sum(!is.na(ga_numeric))
  if (n_ga > 0) {
    cat("Creating PCA plot by gestational age gradient...\n")
    filename <- file.path(dataset_dir, paste0(dataset_id, "_Gestational.Age_gradient_PC1_PC2.svg"))
    tryCatch({
      plot_pca_ga_gradient_single(pca_dataset, pdata_aligned,
                                  axes = c(1, 2), output_file = filename,
                                  output_dir_png = dataset_dir_png,
                                  xlim = xlim, ylim = ylim)
    }, error = function(e) {
      cat("⚠ Error creating GA gradient plot:", conditionMessage(e), "\n")
    })
  } else {
    cat("⚠ No gestational age information - skipping GA gradient plot\n")
  }

  # Plot with sample labels
  cat("Creating PCA plot with sample labels...\n")
  filename <- file.path(dataset_dir, paste0(dataset_id, "_sample_labels_PC1_PC2.svg"))
  tryCatch({
    plot_pca_samples(pca_dataset, pdata_aligned, sample_col,
                    axes = c(1, 2), output_file = filename,
                    output_dir_png = dataset_dir_png,
                    xlim = xlim, ylim = ylim)
  }, error = function(e) {
    cat("⚠ Error creating sample labels plot:", conditionMessage(e), "\n")
  })

  cat("\n✓ Analysis complete for", dataset_id, "\n")
}

cat("\n")
cat("==================================================================\n")
cat("                    ANALYSIS COMPLETE                             \n")
cat("==================================================================\n\n")
cat("All plots saved to:\n")
cat("  SVG:", pca_output_base, "\n")
cat("  PNG:", pca_output_base_png, "\n\n")
