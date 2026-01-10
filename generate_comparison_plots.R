#!/usr/bin/env Rscript

#' Generate logFC-FDR comparison plots for existing test_first_datasets results
#'
#' This script reads the already-generated unfiltered differential expression
#' results and creates comparison plots between baseline and each iteration.
#'
#' @author Expression Integration Pipeline
#' @date 2025-10-22

cat("\n=== Generating logFC-FDR Comparison Plots ===\n\n")

# Load required libraries
library(ggplot2)

#' Plot logFC-FDR comparison between baseline and current iteration
#'
#' @param baseline_difexp Data frame with baseline differential expression results
#' @param current_difexp Data frame with current iteration differential expression results
#' @param output_file Path to save the plot
#' @param title Title for the plot
plot_logfc_fdr_comparison <- function(baseline_difexp, current_difexp, output_file, title = "logFC-FDR Comparison") {

  # Merge baseline and current by gene identifier (using SYMBOL or ENTREZID)
  # Assuming both have SYMBOL column
  gene_col <- if ("SYMBOL" %in% colnames(baseline_difexp)) "SYMBOL" else "ENTREZID"

  baseline_data <- baseline_difexp[, c(gene_col, "logFC", "adj.P.Val")]
  colnames(baseline_data) <- c("gene", "logFC_baseline", "FDR_baseline")

  current_data <- current_difexp[, c(gene_col, "logFC", "adj.P.Val")]
  colnames(current_data) <- c("gene", "logFC_current", "FDR_current")

  # Merge on gene
  merged_data <- merge(baseline_data, current_data, by = "gene", all = FALSE)

  # Convert FDR to -log10(FDR) for better visualization
  merged_data$neglog10FDR_baseline <- -log10(merged_data$FDR_baseline + 1e-300)
  merged_data$neglog10FDR_current <- -log10(merged_data$FDR_current + 1e-300)

  # Cap at reasonable values for visualization
  merged_data$neglog10FDR_baseline <- pmin(merged_data$neglog10FDR_baseline, 50)
  merged_data$neglog10FDR_current <- pmin(merged_data$neglog10FDR_current, 50)

  # Determine point significance
  merged_data$significance <- "Not significant"
  merged_data$significance[merged_data$FDR_baseline < 0.05 | merged_data$FDR_current < 0.05] <- "Significant in one"
  merged_data$significance[merged_data$FDR_baseline < 0.05 & merged_data$FDR_current < 0.05] <- "Significant in both"

  # Calculate change magnitude for coloring lines
  merged_data$logFC_change <- abs(merged_data$logFC_current - merged_data$logFC_baseline)
  merged_data$FDR_change <- abs(merged_data$neglog10FDR_current - merged_data$neglog10FDR_baseline)
  merged_data$total_change <- sqrt(merged_data$logFC_change^2 + merged_data$FDR_change^2)

  # Create the plot
  p <- ggplot(merged_data) +
    # Draw lines connecting baseline and current for each gene
    geom_segment(aes(x = logFC_baseline, y = neglog10FDR_baseline,
                     xend = logFC_current, yend = neglog10FDR_current,
                     color = total_change),
                 alpha = 0.3, linewidth = 0.3) +
    # Baseline points
    geom_point(aes(x = logFC_baseline, y = neglog10FDR_baseline, shape = "Baseline"),
               color = "blue", alpha = 0.5, size = 1.5) +
    # Current iteration points
    geom_point(aes(x = logFC_current, y = neglog10FDR_current, shape = "Current"),
               color = "red", alpha = 0.5, size = 1.5) +
    # Add significance threshold lines
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50", alpha = 0.5) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50", alpha = 0.5) +
    # Color scale for lines
    scale_color_gradient(low = "lightgray", high = "purple", name = "Change\nMagnitude") +
    # Shape scale
    scale_shape_manual(values = c("Baseline" = 16, "Current" = 17), name = "Dataset") +
    # Labels and theme
    labs(
      title = title,
      subtitle = paste0("N genes = ", nrow(merged_data)),
      x = "logFC",
      y = "-log10(FDR)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )

  # Save plot using png device to avoid graphics API issues
  png(output_file, width = 10, height = 8, units = "in", res = 300)
  print(p)
  dev.off()

  cat("  ✓ Saved comparison plot to:", output_file, "\n")

  # Calculate detailed statistics
  # FDR changes
  genes_fdr_increased <- sum(merged_data$FDR_current > merged_data$FDR_baseline)
  genes_fdr_decreased <- sum(merged_data$FDR_current < merged_data$FDR_baseline)

  # Significance changes
  genes_lost_significance <- sum(merged_data$FDR_baseline < 0.05 & merged_data$FDR_current >= 0.05)
  genes_gained_significance <- sum(merged_data$FDR_baseline >= 0.05 & merged_data$FDR_current < 0.05)

  # LogFC changes
  genes_logfc_change_gt1 <- sum(merged_data$logFC_change > 1)

  # LogFC threshold crossings (considering absolute values)
  baseline_above1 <- abs(merged_data$logFC_baseline) > 1
  current_above1 <- abs(merged_data$logFC_current) > 1
  genes_logfc_fell_below1 <- sum(baseline_above1 & !current_above1)
  genes_logfc_rose_above1 <- sum(!baseline_above1 & current_above1)

  # Return summary statistics
  stats <- list(
    n_genes = nrow(merged_data),
    n_sig_baseline = sum(merged_data$FDR_baseline < 0.05),
    n_sig_current = sum(merged_data$FDR_current < 0.05),
    n_sig_both = sum(merged_data$FDR_baseline < 0.05 & merged_data$FDR_current < 0.05),
    mean_logFC_change = mean(merged_data$logFC_change),
    median_logFC_change = median(merged_data$logFC_change),
    genes_fdr_increased = genes_fdr_increased,
    genes_fdr_decreased = genes_fdr_decreased,
    genes_lost_significance = genes_lost_significance,
    genes_gained_significance = genes_gained_significance,
    genes_logfc_change_gt1 = genes_logfc_change_gt1,
    genes_logfc_fell_below1 = genes_logfc_fell_below1,
    genes_logfc_rose_above1 = genes_logfc_rose_above1
  )

  return(stats)
}

# Set paths
test_output_dir <- "output/test_first_datasets"
baseline_file <- file.path(test_output_dir, "BASELINE", "difexp_genes_unfiltered.csv")

# Check if baseline file exists
if (!file.exists(baseline_file)) {
  stop("Baseline file not found: ", baseline_file)
}

# Load baseline data
cat("Loading baseline data...\n")
baseline_difexp <- read.csv(baseline_file, stringsAsFactors = FALSE)
cat("  Baseline genes:", nrow(baseline_difexp), "\n\n")

# Get list of iteration directories (exclude BASELINE)
iteration_dirs <- list.dirs(test_output_dir, full.names = FALSE, recursive = FALSE)
iteration_dirs <- iteration_dirs[iteration_dirs != "BASELINE" & iteration_dirs != ""]

cat("Found", length(iteration_dirs), "iteration(s) to process\n\n")

# Storage for detailed comparison summary
detailed_summary <- data.frame(
  dataset = character(),
  n_genes = integer(),
  n_sig_baseline = integer(),
  n_sig_current = integer(),
  n_sig_both = integer(),
  genes_fdr_increased = integer(),
  genes_fdr_decreased = integer(),
  genes_lost_significance = integer(),
  genes_gained_significance = integer(),
  genes_logfc_change_gt1 = integer(),
  genes_logfc_fell_below1 = integer(),
  genes_logfc_rose_above1 = integer(),
  mean_logFC_change = numeric(),
  median_logFC_change = numeric(),
  stringsAsFactors = FALSE
)

# Process each iteration
for (iter_name in iteration_dirs) {
  cat("Processing:", iter_name, "\n")

  iter_dir <- file.path(test_output_dir, iter_name)
  current_file <- file.path(iter_dir, "difexp_genes_unfiltered.csv")

  # Check if unfiltered file exists
  if (!file.exists(current_file)) {
    cat("  ⚠ Unfiltered file not found, skipping\n\n")
    next
  }

  # Load current iteration data
  current_difexp <- read.csv(current_file, stringsAsFactors = FALSE)
  cat("  Current genes:", nrow(current_difexp), "\n")

  # Generate plot
  plot_file <- file.path(iter_dir, "logFC_FDR_comparison.png")
  plot_title <- paste0("Baseline vs ", iter_name)

  tryCatch({
    plot_stats <- plot_logfc_fdr_comparison(
      baseline_difexp = baseline_difexp,
      current_difexp = current_difexp,
      output_file = plot_file,
      title = plot_title
    )

    cat("  Common genes:", plot_stats$n_genes, "\n")
    cat("  Significant in baseline:", plot_stats$n_sig_baseline, "\n")
    cat("  Significant in current:", plot_stats$n_sig_current, "\n")
    cat("  Significant in both:", plot_stats$n_sig_both, "\n")
    cat("  Genes FDR increased:", plot_stats$genes_fdr_increased, "\n")
    cat("  Genes FDR decreased:", plot_stats$genes_fdr_decreased, "\n")
    cat("  Genes lost significance:", plot_stats$genes_lost_significance, "\n")
    cat("  Genes gained significance:", plot_stats$genes_gained_significance, "\n")
    cat("  Genes |logFC| change > 1:", plot_stats$genes_logfc_change_gt1, "\n")
    cat("  Genes |logFC| fell below 1:", plot_stats$genes_logfc_fell_below1, "\n")
    cat("  Genes |logFC| rose above 1:", plot_stats$genes_logfc_rose_above1, "\n")
    cat("  Mean |logFC| change:", round(plot_stats$mean_logFC_change, 3), "\n")
    cat("  Median |logFC| change:", round(plot_stats$median_logFC_change, 3), "\n")

    # Add to detailed summary
    detailed_summary <- rbind(detailed_summary, data.frame(
      dataset = iter_name,
      n_genes = plot_stats$n_genes,
      n_sig_baseline = plot_stats$n_sig_baseline,
      n_sig_current = plot_stats$n_sig_current,
      n_sig_both = plot_stats$n_sig_both,
      genes_fdr_increased = plot_stats$genes_fdr_increased,
      genes_fdr_decreased = plot_stats$genes_fdr_decreased,
      genes_lost_significance = plot_stats$genes_lost_significance,
      genes_gained_significance = plot_stats$genes_gained_significance,
      genes_logfc_change_gt1 = plot_stats$genes_logfc_change_gt1,
      genes_logfc_fell_below1 = plot_stats$genes_logfc_fell_below1,
      genes_logfc_rose_above1 = plot_stats$genes_logfc_rose_above1,
      mean_logFC_change = plot_stats$mean_logFC_change,
      median_logFC_change = plot_stats$median_logFC_change,
      stringsAsFactors = FALSE
    ))

  }, error = function(e) {
    cat("  ✗ Error generating plot:", conditionMessage(e), "\n")
  })

  cat("\n")
}

# Save detailed summary
summary_file <- file.path(test_output_dir, "detailed_comparison_summary.csv")
write.csv(detailed_summary, summary_file, row.names = FALSE)
cat("\n✓ Saved detailed comparison summary to:", summary_file, "\n")

cat("==================================================================\n")
cat("                    PLOT GENERATION COMPLETE                      \n")
cat("==================================================================\n\n")
