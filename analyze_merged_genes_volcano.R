#!/usr/bin/env Rscript

#' Volcano Plots for Merged vs Excluded Genes
#'
#' For each dataset, create volcano plots showing:
#' - Mean expression (x-axis) vs Variance (y-axis)
#' - Colored by whether genes are merged or excluded
#'
#' This helps visualize the expression-variance relationship
#' and see where excluded genes fall in this space
#'
#' @author Expression Integration Pipeline
#' @date 2025-10-05

cat("\n=== Creating Volcano Plots for Merged vs Excluded Genes ===\n\n")

# Load libraries
library(ggplot2)

# Source pipeline modules to reuse existing functions
source("R/01_data_merging.R")

# Configuration
mapped_dir <- "data/mapped"
pattern <- ".tsv"
output_dir <- "output/gene_volcano_analysis"

# Create output directory
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Get list of expression files
exprs_files <- list.files(mapped_dir, pattern = pattern, full.names = FALSE)

if (length(exprs_files) == 0) {
  stop("No expression files found in ", mapped_dir)
}

cat("Found", length(exprs_files), "expression files\n\n")

# Step 1: Use existing merge function to get common genes
cat("Step 1: Merging all datasets using existing pipeline function...\n")

merged_data <- merge_expression_data(mapped_dir, pattern)
common_genes <- rownames(merged_data)

cat("\nFinal common genes:", length(common_genes), "\n\n")

# Step 2: Create volcano plots for each dataset
cat("Step 2: Creating volcano plots for each dataset...\n\n")

results_summary <- data.frame()

for (file_name in exprs_files) {

  dataset_name <- gsub(".tsv$", "", file_name)
  file_path <- file.path(mapped_dir, file_name)
  cat("Analyzing:", dataset_name, "\n")

  # Read expression data
  exprs <- read.table(file_path, header = TRUE, sep = '\t', row.names = 1)

  # Identify merged and excluded genes
  merged_genes <- intersect(rownames(exprs), common_genes)
  excluded_genes <- setdiff(rownames(exprs), common_genes)

  cat("  Total genes:", nrow(exprs), "\n")
  cat("  Merged genes:", length(merged_genes), "\n")
  cat("  Excluded genes:", length(excluded_genes), "\n")

  # Calculate mean expression and variance per gene
  gene_mean <- rowMeans(exprs, na.rm = TRUE)
  gene_var <- apply(exprs, 1, var, na.rm = TRUE)
  gene_sd <- apply(exprs, 1, sd, na.rm = TRUE)
  gene_cv <- gene_sd / abs(gene_mean)  # Coefficient of variation

  # Create data frame for plotting
  plot_data <- data.frame(
    gene_id = rownames(exprs),
    mean_expression = gene_mean,
    variance = gene_var,
    sd = gene_sd,
    cv = gene_cv,
    category = ifelse(rownames(exprs) %in% merged_genes, "Merged", "Excluded"),
    stringsAsFactors = FALSE
  )

  # Summary stats
  summary_stats <- data.frame(
    dataset = dataset_name,
    total_genes = nrow(exprs),
    merged_genes = length(merged_genes),
    excluded_genes = length(excluded_genes),
    merged_mean_expr = mean(plot_data$mean_expression[plot_data$category == "Merged"], na.rm = TRUE),
    merged_mean_var = mean(plot_data$variance[plot_data$category == "Merged"], na.rm = TRUE),
    excluded_mean_expr = mean(plot_data$mean_expression[plot_data$category == "Excluded"], na.rm = TRUE),
    excluded_mean_var = mean(plot_data$variance[plot_data$category == "Excluded"], na.rm = TRUE)
  )
  results_summary <- rbind(results_summary, summary_stats)

  # Plot 1: Mean vs Variance (classic volcano style)
  cat("  Creating mean vs variance plot...\n")

  # Add small constant to handle zeros for log scale
  plot_data$log_var <- log10(plot_data$variance + 1e-10)

  p <- ggplot(plot_data, aes(x = mean_expression, y = log_var, color = category)) +
    geom_point(alpha = 0.4, size = 1.5) +
    scale_color_manual(values = c("Merged" = "#00BFC4", "Excluded" = "#F8766D")) +
    theme_bw() +
    labs(
      title = paste("Expression-Variance Plot:", dataset_name),
      subtitle = paste("Merged:", length(merged_genes), "genes | Excluded:", length(excluded_genes), "genes"),
      x = "Mean Expression (log2)",
      y = "log10(Variance)",
      color = "Gene Category"
    ) +
    theme(legend.position = "bottom")

  tryCatch({
    ggsave(file.path(output_dir, paste0("volcano_mean_var_", dataset_name, ".png")),
           p, width = 10, height = 8, dpi = 300, device = "png")
  }, error = function(e) {
    png(file.path(output_dir, paste0("volcano_mean_var_", dataset_name, ".png")),
        width = 10*300, height = 8*300, res = 300)
    print(p)
    dev.off()
  })

  # Plot 2: Mean vs Standard Deviation
  cat("  Creating mean vs SD plot...\n")

  p <- ggplot(plot_data, aes(x = mean_expression, y = sd, color = category)) +
    geom_point(alpha = 0.4, size = 1.5) +
    scale_color_manual(values = c("Merged" = "#00BFC4", "Excluded" = "#F8766D")) +
    theme_bw() +
    labs(
      title = paste("Expression-SD Plot:", dataset_name),
      subtitle = paste("Merged:", length(merged_genes), "genes | Excluded:", length(excluded_genes), "genes"),
      x = "Mean Expression (log2)",
      y = "Standard Deviation",
      color = "Gene Category"
    ) +
    theme(legend.position = "bottom")

  tryCatch({
    ggsave(file.path(output_dir, paste0("volcano_mean_sd_", dataset_name, ".png")),
           p, width = 10, height = 8, dpi = 300, device = "png")
  }, error = function(e) {
    png(file.path(output_dir, paste0("volcano_mean_sd_", dataset_name, ".png")),
        width = 10*300, height = 8*300, res = 300)
    print(p)
    dev.off()
  })

  # Plot 3: Mean vs Coefficient of Variation (CV)
  cat("  Creating mean vs CV plot...\n")

  # Filter out infinite CV values
  plot_data_cv <- plot_data[is.finite(plot_data$cv) & plot_data$cv < 100, ]

  if (nrow(plot_data_cv) > 0) {
    p <- ggplot(plot_data_cv, aes(x = mean_expression, y = cv, color = category)) +
      geom_point(alpha = 0.4, size = 1.5) +
      scale_color_manual(values = c("Merged" = "#00BFC4", "Excluded" = "#F8766D")) +
      theme_bw() +
      labs(
        title = paste("Expression-CV Plot:", dataset_name),
        subtitle = paste("Coefficient of Variation (SD/Mean) - Merged:", length(merged_genes), "genes | Excluded:", length(excluded_genes), "genes"),
        x = "Mean Expression (log2)",
        y = "Coefficient of Variation (CV)",
        color = "Gene Category"
      ) +
      theme(legend.position = "bottom")

    tryCatch({
      ggsave(file.path(output_dir, paste0("volcano_mean_cv_", dataset_name, ".png")),
             p, width = 10, height = 8, dpi = 300, device = "png")
    }, error = function(e) {
      png(file.path(output_dir, paste0("volcano_mean_cv_", dataset_name, ".png")),
          width = 10*300, height = 8*300, res = 300)
      print(p)
      dev.off()
    })
  }

  # Plot 4: Density contours overlaid (mean vs variance)
  cat("  Creating density contour plot...\n")

  p <- ggplot(plot_data, aes(x = mean_expression, y = log_var)) +
    geom_point(aes(color = category), alpha = 0.3, size = 1) +
    geom_density_2d(aes(group = category, color = category), alpha = 0.8, linewidth = 0.8) +
    scale_color_manual(values = c("Merged" = "#00BFC4", "Excluded" = "#F8766D")) +
    theme_bw() +
    labs(
      title = paste("Expression-Variance Density:", dataset_name),
      subtitle = "Contour lines show density of genes in each category",
      x = "Mean Expression (log2)",
      y = "log10(Variance)",
      color = "Gene Category"
    ) +
    theme(legend.position = "bottom")

  tryCatch({
    ggsave(file.path(output_dir, paste0("volcano_density_", dataset_name, ".png")),
           p, width = 10, height = 8, dpi = 300, device = "png")
  }, error = function(e) {
    png(file.path(output_dir, paste0("volcano_density_", dataset_name, ".png")),
        width = 10*300, height = 8*300, res = 300)
    print(p)
    dev.off()
  })

  cat("  ✓ Completed plots for", dataset_name, "\n\n")
}

# Save summary statistics
write.csv(results_summary, file.path(output_dir, "volcano_summary.csv"), row.names = FALSE)

# Create overall summary plot - all datasets in one
cat("Step 3: Creating combined summary plot...\n")

# Read all data again for combined plot
all_data <- data.frame()

for (file_name in exprs_files) {
  dataset_name <- gsub(".tsv$", "", file_name)
  file_path <- file.path(mapped_dir, file_name)

  exprs <- read.table(file_path, header = TRUE, sep = '\t', row.names = 1)

  merged_genes <- intersect(rownames(exprs), common_genes)

  gene_mean <- rowMeans(exprs, na.rm = TRUE)
  gene_var <- apply(exprs, 1, var, na.rm = TRUE)

  plot_data <- data.frame(
    dataset = dataset_name,
    gene_id = rownames(exprs),
    mean_expression = gene_mean,
    variance = gene_var,
    log_var = log10(gene_var + 1e-10),
    category = ifelse(rownames(exprs) %in% merged_genes, "Merged", "Excluded"),
    stringsAsFactors = FALSE
  )

  all_data <- rbind(all_data, plot_data)
}

# Combined faceted plot
p <- ggplot(all_data, aes(x = mean_expression, y = log_var, color = category)) +
  geom_point(alpha = 0.2, size = 0.5) +
  facet_wrap(~dataset, ncol = 3) +
  scale_color_manual(values = c("Merged" = "#00BFC4", "Excluded" = "#F8766D")) +
  theme_bw() +
  labs(
    title = "Expression-Variance Plot: All Datasets",
    subtitle = "Comparison of merged vs excluded genes across all datasets",
    x = "Mean Expression (log2)",
    y = "log10(Variance)",
    color = "Gene Category"
  ) +
  theme(legend.position = "bottom")

tryCatch({
  ggsave(file.path(output_dir, "volcano_all_datasets.png"),
         p, width = 16, height = 12, dpi = 300, device = "png")
}, error = function(e) {
  png(file.path(output_dir, "volcano_all_datasets.png"),
      width = 16*300, height = 12*300, res = 300)
  print(p)
  dev.off()
})

cat("\n✓ Analysis complete!\n")
cat("\nResults saved to:", output_dir, "\n")
cat("\nFiles created:\n")
cat("  - volcano_summary.csv (summary statistics)\n")
cat("  - volcano_mean_var_*.png (mean vs variance, one per dataset)\n")
cat("  - volcano_mean_sd_*.png (mean vs SD, one per dataset)\n")
cat("  - volcano_mean_cv_*.png (mean vs CV, one per dataset)\n")
cat("  - volcano_density_*.png (density contours, one per dataset)\n")
cat("  - volcano_all_datasets.png (combined faceted plot)\n\n")

# Print interpretation guide
cat("=== Interpretation Guide ===\n\n")
cat("Mean-Variance Relationship:\n")
cat("  - Low mean, low variance → Lowly expressed, stable genes\n")
cat("  - Low mean, high variance → Lowly expressed, noisy genes (technical noise)\n")
cat("  - High mean, low variance → Highly expressed, stable genes (housekeeping)\n")
cat("  - High mean, high variance → Highly expressed, variable genes (biological signal)\n\n")

cat("If excluded genes cluster in:\n")
cat("  - Low mean region → Filtering out low-expressed genes\n")
cat("  - High variance region → Filtering out noisy genes\n")
cat("  - Random distribution → Annotation/platform differences\n\n")

cat("Coefficient of Variation (CV):\n")
cat("  - Low CV → Stable expression relative to mean\n")
cat("  - High CV → Variable expression relative to mean\n\n")
