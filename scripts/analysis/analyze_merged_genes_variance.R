#!/usr/bin/env Rscript

#' Analyze Expression Variance of Merged vs Excluded Genes
#'
#' For each dataset, compare the expression variance of:
#' - Genes that made it to the merged dataset (common genes)
#' - Genes that were excluded (dataset-specific genes)
#'
#' This helps understand if excluded genes have high/low variance
#'
#' @author Expression Integration Pipeline
#' @date 2025-10-05

cat("\n=== Analyzing Variance of Merged vs Excluded Genes ===\n\n")

# Load libraries
library(ggplot2)

# Source pipeline modules to reuse existing functions
source("R/01_data_merging.R")

# Configuration
mapped_dir <- "data/mapped"
pattern <- ".tsv"
output_dir <- "output/gene_variance_analysis"

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
cat("Step 2: Analyzing variance in each dataset...\n\n")

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

  # Calculate variance per gene (across all samples)
  gene_variance <- apply(exprs, 1, var, na.rm = TRUE)

  # Separate by merged/excluded
  merged_var <- gene_variance[merged_genes]
  excluded_var <- gene_variance[excluded_genes]
  all_var <- gene_variance

  # Calculate statistics
  stats <- data.frame(
    dataset = dataset_name,
    total_genes = nrow(exprs),
    merged_genes = length(merged_genes),
    excluded_genes = length(excluded_genes),

    # Overall statistics
    overall_mean_var = mean(all_var, na.rm = TRUE),
    overall_median_var = median(all_var, na.rm = TRUE),
    overall_sd_var = sd(all_var, na.rm = TRUE),

    # Merged genes statistics
    merged_mean_var = mean(merged_var, na.rm = TRUE),
    merged_median_var = median(merged_var, na.rm = TRUE),
    merged_sd_var = sd(merged_var, na.rm = TRUE),
    merged_min_var = min(merged_var, na.rm = TRUE),
    merged_max_var = max(merged_var, na.rm = TRUE),

    # Excluded genes statistics
    excluded_mean_var = if(length(excluded_var) > 0) mean(excluded_var, na.rm = TRUE) else NA,
    excluded_median_var = if(length(excluded_var) > 0) median(excluded_var, na.rm = TRUE) else NA,
    excluded_sd_var = if(length(excluded_var) > 0) sd(excluded_var, na.rm = TRUE) else NA,
    excluded_min_var = if(length(excluded_var) > 0) min(excluded_var, na.rm = TRUE) else NA,
    excluded_max_var = if(length(excluded_var) > 0) max(excluded_var, na.rm = TRUE) else NA
  )

  results_summary <- rbind(results_summary, stats)

  # Create plotting data
  plot_data <- data.frame(
    dataset = dataset_name,
    gene_id = names(gene_variance),
    variance = gene_variance,
    category = ifelse(names(gene_variance) %in% merged_genes, "Merged", "Excluded"),
    stringsAsFactors = FALSE
  )

  all_plot_data <- rbind(all_plot_data, plot_data)

  cat("  Merged mean variance:", round(stats$merged_mean_var, 4),
      "| Excluded mean variance:", round(stats$excluded_mean_var, 4), "\n")
  cat("  Merged median variance:", round(stats$merged_median_var, 4),
      "| Excluded median variance:", round(stats$excluded_median_var, 4), "\n\n")
}

# Save summary statistics
write.csv(results_summary, file.path(output_dir, "variance_statistics.csv"), row.names = FALSE)
cat("Saved summary statistics to:", file.path(output_dir, "variance_statistics.csv"), "\n\n")

# Step 3: Create visualizations
cat("Step 3: Creating visualizations...\n")

# Plot 1: Density plots for each dataset (log-scale for variance)
cat("  Creating density plots...\n")

for (dataset in unique(all_plot_data$dataset)) {

  subset_data <- all_plot_data[all_plot_data$dataset == dataset, ]

  # Add small constant to handle zeros before log transformation
  subset_data$log_variance <- log10(subset_data$variance + 1e-10)

  # Get statistics for this dataset
  stats <- results_summary[results_summary$dataset == dataset, ]

  p <- ggplot(subset_data, aes(x = log_variance, fill = category)) +
    geom_density(alpha = 0.5) +
    geom_vline(xintercept = log10(stats$overall_mean_var + 1e-10),
               linetype = "solid", color = "black", linewidth = 1) +
    geom_vline(xintercept = log10(stats$overall_median_var + 1e-10),
               linetype = "dashed", color = "black", linewidth = 1) +
    geom_vline(xintercept = log10(stats$merged_mean_var + 1e-10),
               linetype = "solid", color = "#00BFC4", linewidth = 0.8) +
    geom_vline(xintercept = log10(stats$excluded_mean_var + 1e-10),
               linetype = "solid", color = "#F8766D", linewidth = 0.8) +
    scale_fill_manual(values = c("Merged" = "#00BFC4", "Excluded" = "#F8766D")) +
    theme_bw() +
    labs(
      title = paste("Gene Expression Variance Distribution:", dataset),
      subtitle = paste("Merged genes:", stats$merged_genes, "| Excluded genes:", stats$excluded_genes),
      x = "log10(Variance)",
      y = "Density",
      fill = "Gene Category"
    ) +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("Overall mean var: %.4f\nOverall median var: %.4f\nMerged mean var: %.4f\nExcluded mean var: %.4f",
                           stats$overall_mean_var, stats$overall_median_var,
                           stats$merged_mean_var, stats$excluded_mean_var),
             hjust = 1.1, vjust = 1.1, size = 3)

  tryCatch({
    ggsave(file.path(output_dir, paste0("density_", dataset, ".png")),
           p, width = 10, height = 6, dpi = 300, device = "png")
  }, error = function(e) {
    png(file.path(output_dir, paste0("density_", dataset, ".png")),
        width = 10*300, height = 6*300, res = 300)
    print(p)
    dev.off()
  })
}

# Plot 2: Boxplots comparing merged vs excluded across all datasets (log-scale)
cat("  Creating comparison boxplots...\n")

all_plot_data$log_variance <- log10(all_plot_data$variance + 1e-10)

p <- ggplot(all_plot_data, aes(x = dataset, y = log_variance, fill = category)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3) +
  scale_fill_manual(values = c("Merged" = "#00BFC4", "Excluded" = "#F8766D")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Expression Variance: Merged vs Excluded Genes Across Datasets",
    x = "Dataset",
    y = "log10(Variance)",
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

# Plot 3: Summary plot - mean variance comparison
cat("  Creating summary plot...\n")

summary_long <- data.frame(
  dataset = rep(results_summary$dataset, 3),
  category = rep(c("Merged", "Excluded", "Overall"), each = nrow(results_summary)),
  mean_variance = c(results_summary$merged_mean_var, results_summary$excluded_mean_var, results_summary$overall_mean_var)
)

# Log transform for better visualization
summary_long$log_mean_variance <- log10(summary_long$mean_variance + 1e-10)

p <- ggplot(summary_long, aes(x = dataset, y = log_mean_variance, color = category, group = category)) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("Merged" = "#00BFC4", "Excluded" = "#F8766D", "Overall" = "black")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Mean Variance: Merged vs Excluded vs Overall",
    subtitle = "Comparison of mean variance across all datasets",
    x = "Dataset",
    y = "log10(Mean Variance)",
    color = "Gene Category"
  )

tryCatch({
  ggsave(file.path(output_dir, "mean_variance_comparison.png"),
         p, width = 12, height = 6, dpi = 300, device = "png")
}, error = function(e) {
  png(file.path(output_dir, "mean_variance_comparison.png"),
      width = 12*300, height = 6*300, res = 300)
  print(p)
  dev.off()
})

# Plot 4: Median variance comparison
summary_long_median <- data.frame(
  dataset = rep(results_summary$dataset, 3),
  category = rep(c("Merged", "Excluded", "Overall"), each = nrow(results_summary)),
  median_variance = c(results_summary$merged_median_var, results_summary$excluded_median_var, results_summary$overall_median_var)
)

summary_long_median$log_median_variance <- log10(summary_long_median$median_variance + 1e-10)

p <- ggplot(summary_long_median, aes(x = dataset, y = log_median_variance, color = category, group = category)) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("Merged" = "#00BFC4", "Excluded" = "#F8766D", "Overall" = "black")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Median Variance: Merged vs Excluded vs Overall",
    subtitle = "Comparison of median variance across all datasets",
    x = "Dataset",
    y = "log10(Median Variance)",
    color = "Gene Category"
  )

tryCatch({
  ggsave(file.path(output_dir, "median_variance_comparison.png"),
         p, width = 12, height = 6, dpi = 300, device = "png")
}, error = function(e) {
  png(file.path(output_dir, "median_variance_comparison.png"),
      width = 12*300, height = 6*300, res = 300)
  print(p)
  dev.off()
})

# Plot 5: Difference plot (merged - excluded variance)
results_summary$variance_diff <- results_summary$merged_mean_var - results_summary$excluded_mean_var

p <- ggplot(results_summary, aes(x = dataset, y = variance_diff)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Variance Difference: Merged - Excluded Genes",
    subtitle = "Positive values = merged genes have higher variance",
    x = "Dataset",
    y = "Difference in Mean Variance"
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
cat("  - variance_statistics.csv (summary table)\n")
cat("  - density_*.png (one per dataset)\n")
cat("  - boxplot_comparison.png\n")
cat("  - mean_variance_comparison.png\n")
cat("  - median_variance_comparison.png\n")
cat("  - difference_plot.png\n\n")

# Print interpretation guide
cat("=== Interpretation Guide ===\n\n")
cat("If excluded genes have:\n")
cat("  - LOWER variance than merged genes:\n")
cat("    → Excluded genes are stable/constant across samples\n")
cat("    → Less informative for differential expression\n\n")
cat("  - HIGHER variance than merged genes:\n")
cat("    → Excluded genes are noisy/variable\n")
cat("    → May be dataset-specific technical variation\n\n")
cat("  - SIMILAR variance to merged genes:\n")
cat("    → Variance not driving exclusion\n")
cat("    → Genes excluded due to annotation differences\n\n")
