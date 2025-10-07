#!/usr/bin/env Rscript

#' Analyze Expression Levels of Merged vs Excluded Genes
#'
#' For each dataset, compare the expression distribution of:
#' - Genes that made it to the merged dataset (common genes)
#' - Genes that were excluded (dataset-specific genes)
#'
#' This helps understand if excluded genes are low/high expressed
#'
#' @author Expression Integration Pipeline
#' @date 2025-10-05

cat("\n=== Analyzing Merged vs Excluded Genes ===\n\n")

# Load libraries
library(ggplot2)

# Source pipeline modules to reuse existing functions
source("R/01_data_merging.R")

# Configuration
mapped_dir <- "data/mapped"
pattern <- ".tsv"
output_dir <- "output/gene_inclusion_analysis"

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

# Step 2: Analyze each dataset
cat("Step 2: Analyzing each dataset...\n\n")

results_summary <- data.frame()
all_plot_data <- data.frame()

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

  # Calculate average expression per gene (across all samples)
  gene_avg_exprs <- rowMeans(exprs, na.rm = TRUE)

  # Separate by merged/excluded
  merged_avg <- gene_avg_exprs[merged_genes]
  excluded_avg <- gene_avg_exprs[excluded_genes]
  all_avg <- gene_avg_exprs

  # Calculate statistics
  stats <- data.frame(
    dataset = dataset_name,
    total_genes = nrow(exprs),
    merged_genes = length(merged_genes),
    excluded_genes = length(excluded_genes),

    # Overall statistics
    overall_mean = mean(all_avg, na.rm = TRUE),
    overall_median = median(all_avg, na.rm = TRUE),

    # Merged genes statistics
    merged_mean = mean(merged_avg, na.rm = TRUE),
    merged_median = median(merged_avg, na.rm = TRUE),
    merged_min = min(merged_avg, na.rm = TRUE),
    merged_max = max(merged_avg, na.rm = TRUE),

    # Excluded genes statistics
    excluded_mean = if(length(excluded_avg) > 0) mean(excluded_avg, na.rm = TRUE) else NA,
    excluded_median = if(length(excluded_avg) > 0) median(excluded_avg, na.rm = TRUE) else NA,
    excluded_min = if(length(excluded_avg) > 0) min(excluded_avg, na.rm = TRUE) else NA,
    excluded_max = if(length(excluded_avg) > 0) max(excluded_avg, na.rm = TRUE) else NA
  )

  results_summary <- rbind(results_summary, stats)

  # Create plotting data
  plot_data <- data.frame(
    dataset = dataset_name,
    gene_id = names(gene_avg_exprs),
    avg_expression = gene_avg_exprs,
    category = ifelse(names(gene_avg_exprs) %in% merged_genes, "Merged", "Excluded"),
    stringsAsFactors = FALSE
  )

  all_plot_data <- rbind(all_plot_data, plot_data)

  cat("  Merged mean:", round(stats$merged_mean, 2),
      "| Excluded mean:", round(stats$excluded_mean, 2), "\n")
  cat("  Merged median:", round(stats$merged_median, 2),
      "| Excluded median:", round(stats$excluded_median, 2), "\n\n")
}

# Save summary statistics
write.csv(results_summary, file.path(output_dir, "merge_statistics.csv"), row.names = FALSE)
cat("Saved summary statistics to:", file.path(output_dir, "merge_statistics.csv"), "\n\n")

# Step 3: Create visualizations
cat("Step 3: Creating visualizations...\n")

# Plot 1: Density plots for each dataset
cat("  Creating density plots...\n")

for (dataset in unique(all_plot_data$dataset)) {

  subset_data <- all_plot_data[all_plot_data$dataset == dataset, ]

  # Get statistics for this dataset
  stats <- results_summary[results_summary$dataset == dataset, ]

  p <- ggplot(subset_data, aes(x = avg_expression, fill = category)) +
    geom_density(alpha = 0.5) +
    geom_vline(xintercept = stats$overall_mean, linetype = "solid", color = "black", linewidth = 1) +
    geom_vline(xintercept = stats$overall_median, linetype = "dashed", color = "black", linewidth = 1) +
    geom_vline(xintercept = stats$merged_mean, linetype = "solid", color = "#00BFC4", linewidth = 0.8) +
    geom_vline(xintercept = stats$excluded_mean, linetype = "solid", color = "#F8766D", linewidth = 0.8) +
    scale_fill_manual(values = c("Merged" = "#00BFC4", "Excluded" = "#F8766D")) +
    theme_bw() +
    labs(
      title = paste("Gene Expression Distribution:", dataset),
      subtitle = paste("Merged genes:", stats$merged_genes, "| Excluded genes:", stats$excluded_genes),
      x = "Average Expression (log2)",
      y = "Density",
      fill = "Gene Category"
    ) +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("Overall mean: %.2f\nOverall median: %.2f\nMerged mean: %.2f\nExcluded mean: %.2f",
                           stats$overall_mean, stats$overall_median,
                           stats$merged_mean, stats$excluded_mean),
             hjust = 1.1, vjust = 1.1, size = 3)

  tryCatch({
    ggsave(file.path(output_dir, paste0("density_", dataset, ".png")),
           p, width = 10, height = 6, dpi = 300, device = "png")
  }, error = function(e) {
    # Fallback: use basic png device
    png(file.path(output_dir, paste0("density_", dataset, ".png")),
        width = 10*300, height = 6*300, res = 300)
    print(p)
    dev.off()
  })
}

# Plot 2: Boxplots comparing merged vs excluded across all datasets
cat("  Creating comparison boxplots...\n")

p <- ggplot(all_plot_data, aes(x = dataset, y = avg_expression, fill = category)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3) +
  scale_fill_manual(values = c("Merged" = "#00BFC4", "Excluded" = "#F8766D")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Expression Levels: Merged vs Excluded Genes Across Datasets",
    x = "Dataset",
    y = "Average Expression (log2)",
    fill = "Gene Category"
  )

tryCatch({
  ggsave(file.path(output_dir, "boxplot_comparison.png"),
         p, width = 12, height = 6, dpi = 300, device = "png")
}, error = function(e) {
  png(file.path(output_dir, "boxplot_comparison.png"),
      width = 12*300, height = 6*300, res = 300)
  print(p)
  dev.off()
})

# Plot 3: Summary plot - mean expression comparison
cat("  Creating summary plot...\n")

# Reshape for plotting - include all three: Merged, Excluded, and Overall
summary_long <- data.frame(
  dataset = rep(results_summary$dataset, 3),
  category = rep(c("Merged", "Excluded", "Overall"), each = nrow(results_summary)),
  mean_expression = c(results_summary$merged_mean, results_summary$excluded_mean, results_summary$overall_mean)
)

p <- ggplot(summary_long, aes(x = dataset, y = mean_expression, color = category, group = category)) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("Merged" = "#00BFC4", "Excluded" = "#F8766D", "Overall" = "black")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Mean Expression: Merged vs Excluded vs Overall",
    subtitle = "Comparison of mean expression across all datasets",
    x = "Dataset",
    y = "Mean Expression (log2)",
    color = "Gene Category"
  )

tryCatch({
  ggsave(file.path(output_dir, "mean_comparison.png"),
         p, width = 12, height = 6, dpi = 300, device = "png")
}, error = function(e) {
  png(file.path(output_dir, "mean_comparison.png"),
      width = 12*300, height = 6*300, res = 300)
  print(p)
  dev.off()
})

# Plot 4: Median expression comparison
summary_long_median <- data.frame(
  dataset = rep(results_summary$dataset, 3),
  category = rep(c("Merged", "Excluded", "Overall"), each = nrow(results_summary)),
  median_expression = c(results_summary$merged_median, results_summary$excluded_median, results_summary$overall_median)
)

p <- ggplot(summary_long_median, aes(x = dataset, y = median_expression, color = category, group = category)) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("Merged" = "#00BFC4", "Excluded" = "#F8766D", "Overall" = "black")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Median Expression: Merged vs Excluded vs Overall",
    subtitle = "Comparison of median expression across all datasets",
    x = "Dataset",
    y = "Median Expression (log2)",
    color = "Gene Category"
  )

tryCatch({
  ggsave(file.path(output_dir, "median_comparison.png"),
         p, width = 12, height = 6, dpi = 300, device = "png")
}, error = function(e) {
  png(file.path(output_dir, "median_comparison.png"),
      width = 12*300, height = 6*300, res = 300)
  print(p)
  dev.off()
})

# Plot 5: Difference plot (merged - excluded)
results_summary$expression_diff <- results_summary$merged_mean - results_summary$excluded_mean

p <- ggplot(results_summary, aes(x = dataset, y = expression_diff)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Expression Difference: Merged - Excluded Genes",
    subtitle = "Positive values = merged genes have higher expression",
    x = "Dataset",
    y = "Difference in Mean Expression (log2)"
  )

tryCatch({
  ggsave(file.path(output_dir, "difference_plot.png"),
         p, width = 12, height = 6, dpi = 300, device = "png")
}, error = function(e) {
  png(file.path(output_dir, "difference_plot.png"),
      width = 12*300, height = 6*300, res = 300)
  print(p)
  dev.off()
})

cat("\n✓ Analysis complete!\n")
cat("\nResults saved to:", output_dir, "\n")
cat("\nFiles created:\n")
cat("  - merge_statistics.csv (summary table)\n")
cat("  - density_*.png (one per dataset)\n")
cat("  - boxplot_comparison.png\n")
cat("  - mean_comparison.png\n")
cat("  - median_comparison.png\n")
cat("  - difference_plot.png\n\n")

# Print interpretation guide
cat("=== Interpretation Guide ===\n\n")
cat("If excluded genes have:\n")
cat("  - LOWER expression than merged genes:\n")
cat("    → Low-expressed genes are dataset-specific (less reliable)\n")
cat("    → Merging filters out noise\n\n")
cat("  - HIGHER expression than merged genes:\n")
cat("    → Dataset has unique highly-expressed genes\n")
cat("    → May indicate dataset-specific biology or technical issues\n\n")
cat("  - SIMILAR expression to merged genes:\n")
cat("    → Exclusion is random (not expression-driven)\n")
cat("    → Genes excluded due to annotation differences\n\n")
