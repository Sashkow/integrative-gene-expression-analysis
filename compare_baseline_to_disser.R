#!/usr/bin/env Rscript

#' Compare Baseline Results to Dissertation Data
#'
#' This script compares the baseline differential expression results
#' from test_first_datasets with the dissertation data.
#'
#' @author Expression Integration Pipeline
#' @date 2025-10-26

cat("\n=== Comparing Baseline to Dissertation Data ===\n\n")

# Load required libraries
library(readxl)
library(ggplot2)
library(VennDiagram)

# Define paths
baseline_dir <- "output/test_first_datasets/BASELINE"
disser_file <- "data/disser/lykhenko_Supplement1_sorted_include_present.xlsx"
output_dir <- "output/baseline_vs_disser_comparison"

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
cat("Created output directory:", output_dir, "\n\n")

# ===================================================================
# STEP 1: LOAD BASELINE DATA
# ===================================================================

cat("Loading baseline data...\n")

# Load baseline differential expression results
baseline_difexp_file <- file.path(baseline_dir, "difexp_genes.csv")
baseline_difexp_all_file <- file.path(baseline_dir, "difexp_genes_unfiltered.csv")

if (!file.exists(baseline_difexp_file)) {
  stop("Baseline differential expression file not found: ", baseline_difexp_file,
       "\nPlease run test_first_datasets.R first to generate baseline results.")
}

baseline_difexp <- read.csv(baseline_difexp_file, stringsAsFactors = FALSE)
baseline_difexp_all <- read.csv(baseline_difexp_all_file, stringsAsFactors = FALSE)

cat("✓ Loaded baseline results:\n")
cat("  Filtered DEGs (|logFC| > 1, p < 0.05):", nrow(baseline_difexp), "\n")
cat("  All genes:", nrow(baseline_difexp_all), "\n\n")

# ===================================================================
# STEP 2: LOAD DISSERTATION DATA
# ===================================================================

cat("Loading dissertation data...\n")

if (!file.exists(disser_file)) {
  stop("Dissertation file not found: ", disser_file)
}

# Read first sheet
disser_data <- read_excel(disser_file, sheet = 1)

cat("✓ Loaded dissertation data:\n")
cat("  Total genes:", nrow(disser_data), "\n")
cat("  Columns:", paste(colnames(disser_data), collapse = ", "), "\n\n")

# ===================================================================
# STEP 3: PREPARE DATA FOR COMPARISON
# ===================================================================

cat("Preparing data for comparison...\n")

# Identify common gene identifier column
# Check what columns are available in both datasets
cat("Baseline columns:", paste(colnames(baseline_difexp), collapse = ", "), "\n")
cat("Dissertation columns:", paste(colnames(disser_data), collapse = ", "), "\n\n")

# Use ENTREZID or SYMBOL for comparison (adjust based on actual column names)
gene_col_baseline <- if ("ENTREZID" %in% colnames(baseline_difexp)) "ENTREZID" else "SYMBOL"
gene_col_disser <- if ("ENTREZID" %in% colnames(disser_data)) "ENTREZID" else
                   if ("SYMBOL" %in% colnames(disser_data)) "SYMBOL" else colnames(disser_data)[1]

cat("Using comparison column:\n")
cat("  Baseline:", gene_col_baseline, "\n")
cat("  Dissertation:", gene_col_disser, "\n\n")

# Extract gene lists
# Use ALL baseline genes (unfiltered) for comparison
baseline_genes <- baseline_difexp_all[[gene_col_baseline]]
baseline_genes_filtered <- baseline_difexp[[gene_col_baseline]]
disser_genes <- disser_data[[gene_col_disser]]

# Remove NA values
baseline_genes <- baseline_genes[!is.na(baseline_genes)]
baseline_genes_filtered <- baseline_genes_filtered[!is.na(baseline_genes_filtered)]
disser_genes <- disser_genes[!is.na(disser_genes)]

cat("Gene counts after removing NAs:\n")
cat("  Baseline (all genes):", length(baseline_genes), "\n")
cat("  Baseline (filtered, |logFC| > 1, adj.P.Val < 0.05):", length(baseline_genes_filtered), "\n")
cat("  Dissertation:", length(disser_genes), "\n\n")

# ===================================================================
# STEP 4: VENN DIAGRAM COMPARISON
# ===================================================================

cat("Creating Venn diagram...\n")

# Find overlaps
genes_common <- intersect(baseline_genes, disser_genes)
genes_baseline_only <- setdiff(baseline_genes, disser_genes)
genes_disser_only <- setdiff(disser_genes, baseline_genes)

cat("Overlap statistics:\n")
cat("  Common genes:", length(genes_common), "\n")
cat("  Baseline only:", length(genes_baseline_only), "\n")
cat("  Dissertation only:", length(genes_disser_only), "\n\n")

# Calculate percentages
pct_overlap_baseline <- round(100 * length(genes_common) / length(baseline_genes), 1)
pct_overlap_disser <- round(100 * length(genes_common) / length(disser_genes), 1)

cat("Overlap percentages:\n")
cat("  Of baseline genes:", pct_overlap_baseline, "%\n")
cat("  Of dissertation genes:", pct_overlap_disser, "%\n\n")

# Create Venn diagram
venn_plot <- venn.diagram(
  x = list(
    Baseline = baseline_genes,
    Dissertation = disser_genes
  ),
  filename = NULL,
  fill = c("#440154ff", "#21908dff"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.dist = c(0.055, 0.055)
)

# Save Venn diagram
venn_file <- file.path(output_dir, "venn_diagram.png")
png(venn_file, width = 8, height = 8, units = "in", res = 300)
grid::grid.draw(venn_plot)
dev.off()

cat("✓ Saved Venn diagram to:", venn_file, "\n\n")

# ===================================================================
# STEP 5: COMPARE LOGFC VALUES FOR COMMON GENES
# ===================================================================

cat("Comparing logFC values for common genes...\n")

# Check if both datasets have logFC
if ("logFC" %in% colnames(baseline_difexp_all) && "logFC" %in% colnames(disser_data)) {

  # Check if STRING_id is available in both datasets for composite key merging
  use_string_id <- "STRING_id" %in% colnames(baseline_difexp_all) &&
                   "STRING_id" %in% colnames(disser_data)

  if (use_string_id) {
    cat("  Using composite key (gene + STRING_id) for merging to handle multiple isoforms\n")

    # Merge data for common genes using composite key
    baseline_for_merge <- baseline_difexp_all[, c(gene_col_baseline, "STRING_id", "logFC", "adj.P.Val")]
    colnames(baseline_for_merge) <- c("gene", "STRING_id", "logFC_baseline", "FDR_baseline")

    disser_for_merge <- disser_data[, c(gene_col_disser, "STRING_id", "logFC")]
    colnames(disser_for_merge) <- c("gene", "STRING_id", "logFC_disser")
  } else {
    cat("  STRING_id not available in both datasets, deduplicating by gene ID\n")

    # Deduplicate by keeping only the first occurrence of each gene
    baseline_for_merge <- baseline_difexp_all[, c(gene_col_baseline, "logFC", "adj.P.Val")]
    colnames(baseline_for_merge) <- c("gene", "logFC_baseline", "FDR_baseline")
    baseline_for_merge <- baseline_for_merge[!duplicated(baseline_for_merge$gene), ]

    disser_for_merge <- disser_data[, c(gene_col_disser, "logFC")]
    colnames(disser_for_merge) <- c("gene", "logFC_disser")
    disser_for_merge <- disser_for_merge[!duplicated(disser_for_merge$gene), ]

    cat("  After deduplication: baseline =", nrow(baseline_for_merge),
        "genes, dissertation =", nrow(disser_for_merge), "genes\n")
  }

  # Also get FDR from dissertation if available
  if (use_string_id) {
    # When using STRING_id, extract FDR with STRING_id
    if ("adj.P.Val" %in% colnames(disser_data)) {
      idx <- match(paste(disser_for_merge$gene, disser_for_merge$STRING_id),
                   paste(disser_data[[gene_col_disser]], disser_data$STRING_id))
      disser_for_merge$FDR_disser <- as.numeric(disser_data$adj.P.Val[idx])
    } else if ("P.Value" %in% colnames(disser_data)) {
      idx <- match(paste(disser_for_merge$gene, disser_for_merge$STRING_id),
                   paste(disser_data[[gene_col_disser]], disser_data$STRING_id))
      disser_for_merge$FDR_disser <- as.numeric(disser_data$P.Value[idx])
    } else if ("FDR" %in% colnames(disser_data)) {
      idx <- match(paste(disser_for_merge$gene, disser_for_merge$STRING_id),
                   paste(disser_data[[gene_col_disser]], disser_data$STRING_id))
      disser_for_merge$FDR_disser <- as.numeric(disser_data$FDR[idx])
    }
    # Merge using composite key
    merged_data <- merge(baseline_for_merge, disser_for_merge, by = c("gene", "STRING_id"), all = FALSE)
  } else {
    # When not using STRING_id, FDR already extracted during deduplication
    if ("adj.P.Val" %in% colnames(disser_data)) {
      idx <- match(disser_for_merge$gene, disser_data[[gene_col_disser]])
      disser_for_merge$FDR_disser <- as.numeric(disser_data$adj.P.Val[idx])
    } else if ("P.Value" %in% colnames(disser_data)) {
      idx <- match(disser_for_merge$gene, disser_data[[gene_col_disser]])
      disser_for_merge$FDR_disser <- as.numeric(disser_data$P.Value[idx])
    } else if ("FDR" %in% colnames(disser_data)) {
      idx <- match(disser_for_merge$gene, disser_data[[gene_col_disser]])
      disser_for_merge$FDR_disser <- as.numeric(disser_data$FDR[idx])
    }
    # Merge using single key
    merged_data <- merge(baseline_for_merge, disser_for_merge, by = "gene", all = FALSE)
  }

  cat("  Merged genes for comparison:", nrow(merged_data), "\n")

  # Calculate correlation
  cor_pearson <- cor(merged_data$logFC_baseline, merged_data$logFC_disser, use = "complete.obs")
  cor_spearman <- cor(merged_data$logFC_baseline, merged_data$logFC_disser,
                      method = "spearman", use = "complete.obs")

  cat("  Pearson correlation:", round(cor_pearson, 3), "\n")
  cat("  Spearman correlation:", round(cor_spearman, 3), "\n\n")

  # Create scatter plot
  p <- ggplot(merged_data, aes(x = logFC_disser, y = logFC_baseline)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    geom_smooth(method = "lm", color = "blue", se = TRUE) +
    labs(
      title = "logFC Comparison: Baseline vs Dissertation",
      subtitle = paste0("N = ", nrow(merged_data),
                       " | Pearson r = ", round(cor_pearson, 3),
                       " | Spearman r = ", round(cor_spearman, 3)),
      x = "logFC (Dissertation)",
      y = "logFC (Baseline)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    )

  scatter_file <- file.path(output_dir, "logFC_scatter.png")
  png(scatter_file, width = 8, height = 8, units = "in", res = 300)
  print(p)
  dev.off()

  cat("✓ Saved scatter plot to:", scatter_file, "\n\n")

  # Save merged comparison data
  merged_file <- file.path(output_dir, "common_genes_comparison.csv")
  write.csv(merged_data, merged_file, row.names = FALSE)
  cat("✓ Saved common genes comparison to:", merged_file, "\n\n")

  # Create volcano-style logFC-FDR comparison plot
  if ("FDR_disser" %in% colnames(merged_data) && "FDR_baseline" %in% colnames(merged_data)) {
    cat("Creating logFC-FDR comparison plot...\n")

    # Convert FDR to -log10(FDR) for better visualization
    merged_data$neglog10FDR_baseline <- -log10(merged_data$FDR_baseline + 1e-300)
    merged_data$neglog10FDR_disser <- -log10(merged_data$FDR_disser + 1e-300)

    # Cap at reasonable values for visualization
    merged_data$neglog10FDR_baseline <- pmin(merged_data$neglog10FDR_baseline, 50)
    merged_data$neglog10FDR_disser <- pmin(merged_data$neglog10FDR_disser, 50)

    # Determine point significance
    merged_data$significance <- "Not significant"
    merged_data$significance[merged_data$FDR_baseline < 0.05 | merged_data$FDR_disser < 0.05] <- "Significant in one"
    merged_data$significance[merged_data$FDR_baseline < 0.05 & merged_data$FDR_disser < 0.05] <- "Significant in both"

    # Calculate change magnitude for coloring lines
    merged_data$logFC_change <- abs(merged_data$logFC_baseline - merged_data$logFC_disser)
    merged_data$FDR_change <- abs(merged_data$neglog10FDR_baseline - merged_data$neglog10FDR_disser)
    merged_data$total_change <- sqrt(merged_data$logFC_change^2 + merged_data$FDR_change^2)

    # Create the plot
    p_volcano <- ggplot(merged_data) +
      # Draw lines connecting baseline and dissertation for each gene
      geom_segment(aes(x = logFC_baseline, y = neglog10FDR_baseline,
                       xend = logFC_disser, yend = neglog10FDR_disser,
                       color = total_change),
                   alpha = 0.3, linewidth = 0.3) +
      # Baseline points
      geom_point(aes(x = logFC_baseline, y = neglog10FDR_baseline, shape = "Baseline"),
                 color = "blue", alpha = 0.5, size = 1.5) +
      # Dissertation points
      geom_point(aes(x = logFC_disser, y = neglog10FDR_disser, shape = "Dissertation"),
                 color = "red", alpha = 0.5, size = 1.5) +
      # Add significance threshold lines
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50", alpha = 0.5) +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50", alpha = 0.5) +
      # Color scale for lines
      scale_color_gradient(low = "lightgray", high = "purple", name = "Change\nMagnitude") +
      # Shape scale
      scale_shape_manual(values = c("Baseline" = 16, "Dissertation" = 17), name = "Dataset") +
      # Labels and theme
      labs(
        title = "logFC-FDR Comparison: Baseline vs Dissertation",
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

    # Save plot
    volcano_file <- file.path(output_dir, "logFC_FDR_comparison.png")
    png(volcano_file, width = 10, height = 8, units = "in", res = 300)
    print(p_volcano)
    dev.off()

    cat("✓ Saved logFC-FDR comparison plot to:", volcano_file, "\n\n")

    # Calculate detailed statistics
    # FDR changes
    genes_fdr_increased <- sum(merged_data$FDR_baseline > merged_data$FDR_disser, na.rm = TRUE)
    genes_fdr_decreased <- sum(merged_data$FDR_baseline < merged_data$FDR_disser, na.rm = TRUE)

    # Significance changes
    genes_lost_significance <- sum(merged_data$FDR_baseline < 0.05 & merged_data$FDR_disser >= 0.05, na.rm = TRUE)
    genes_gained_significance <- sum(merged_data$FDR_baseline >= 0.05 & merged_data$FDR_disser < 0.05, na.rm = TRUE)

    # LogFC changes
    genes_logfc_change_gt1 <- sum(merged_data$logFC_change > 1, na.rm = TRUE)

    # LogFC threshold crossings (considering absolute values)
    baseline_above1 <- abs(merged_data$logFC_baseline) > 1
    disser_above1 <- abs(merged_data$logFC_disser) > 1
    genes_logfc_fell_below1 <- sum(baseline_above1 & !disser_above1, na.rm = TRUE)
    genes_logfc_rose_above1 <- sum(!baseline_above1 & disser_above1, na.rm = TRUE)

    cat("LogFC-FDR comparison statistics:\n")
    cat("  Genes with FDR increase (baseline -> disser):", genes_fdr_increased, "\n")
    cat("  Genes with FDR decrease (baseline -> disser):", genes_fdr_decreased, "\n")
    cat("  Genes losing significance:", genes_lost_significance, "\n")
    cat("  Genes gaining significance:", genes_gained_significance, "\n")
    cat("  Genes with |logFC| change > 1:", genes_logfc_change_gt1, "\n")
    cat("  Genes with |logFC| falling below 1:", genes_logfc_fell_below1, "\n")
    cat("  Genes with |logFC| rising above 1:", genes_logfc_rose_above1, "\n\n")
  }

} else {
  cat("⚠ Cannot compare logFC values - column not found in one or both datasets\n\n")
}

# ===================================================================
# STEP 6: DIRECTION OF CHANGE COMPARISON
# ===================================================================

cat("Analyzing direction of change for common genes...\n")

if (exists("merged_data") && nrow(merged_data) > 0) {

  # Classify direction
  merged_data$direction_baseline <- ifelse(merged_data$logFC_baseline > 0, "Up", "Down")
  merged_data$direction_disser <- ifelse(merged_data$logFC_disser > 0, "Up", "Down")
  merged_data$direction_match <- merged_data$direction_baseline == merged_data$direction_disser

  # Count matches
  n_match <- sum(merged_data$direction_match, na.rm = TRUE)
  n_total <- nrow(merged_data)
  pct_match <- round(100 * n_match / n_total, 1)

  cat("Direction concordance:\n")
  cat("  Same direction:", n_match, "(", pct_match, "%)\n")
  cat("  Opposite direction:", n_total - n_match, "(", 100 - pct_match, "%)\n\n")

  # Create contingency table
  direction_table <- table(merged_data$direction_baseline, merged_data$direction_disser)
  cat("Contingency table:\n")
  print(direction_table)
  cat("\n")

}

# ===================================================================
# STEP 7: CREATE DETAILED COMPARISON TABLE
# ===================================================================

cat("Creating detailed comparison table...\n")

if (exists("merged_data") && nrow(merged_data) > 0) {

  # Find genes in dissertation but NOT in baseline at all
  genes_in_disser_not_in_baseline <- setdiff(disser_genes, baseline_genes)

  # Find genes that are significant in baseline (|logFC| > 1, adj.P.Val < 0.05) but NOT in dissertation
  genes_baseline_sig_not_in_disser <- setdiff(baseline_genes_filtered, disser_genes)

  cat("  Genes in dissertation but NOT in baseline:", length(genes_in_disser_not_in_baseline), "\n")
  cat("  Genes significant in baseline (|logFC| > 1, adj.P.Val < 0.05) but NOT in dissertation:",
      length(genes_baseline_sig_not_in_disser), "\n\n")

  # Create comprehensive comparison table with all genes
  comparison_table <- merged_data

  # Add classification column
  comparison_table$classification <- "Common to both"

  # Mark genes that are significant in baseline
  comparison_table$baseline_significant <-
    (abs(comparison_table$logFC_baseline) > 1) &
    (comparison_table$FDR_baseline < 0.05)

  # Mark genes that are significant in dissertation (if FDR available)
  if ("FDR_disser" %in% colnames(comparison_table)) {
    comparison_table$disser_significant <-
      (abs(comparison_table$logFC_disser) > 1) &
      (comparison_table$FDR_disser < 0.05)
  }

  # Check if STRING_id is present in comparison_table
  has_string_id <- "STRING_id" %in% colnames(comparison_table)

  # Add genes in dissertation but not in baseline
  if (length(genes_in_disser_not_in_baseline) > 0) {
    disser_for_merge_extra <- disser_data[disser_data[[gene_col_disser]] %in% genes_in_disser_not_in_baseline, ]

    # Create rows for these genes
    if (has_string_id) {
      extra_rows <- data.frame(
        gene = disser_for_merge_extra[[gene_col_disser]],
        STRING_id = disser_for_merge_extra$STRING_id,
        logFC_baseline = NA,
        FDR_baseline = NA,
        logFC_disser = disser_for_merge_extra$logFC,
        stringsAsFactors = FALSE
      )
    } else {
      extra_rows <- data.frame(
        gene = disser_for_merge_extra[[gene_col_disser]],
        logFC_baseline = NA,
        FDR_baseline = NA,
        logFC_disser = disser_for_merge_extra$logFC,
        stringsAsFactors = FALSE
      )
    }

    # Add FDR_disser if available
    if ("adj.P.Val" %in% colnames(disser_for_merge_extra)) {
      extra_rows$FDR_disser <- disser_for_merge_extra$adj.P.Val
    } else if ("P.Value" %in% colnames(disser_for_merge_extra)) {
      extra_rows$FDR_disser <- disser_for_merge_extra$P.Value
    } else if ("FDR" %in% colnames(disser_for_merge_extra)) {
      extra_rows$FDR_disser <- disser_for_merge_extra$FDR
    }

    extra_rows$classification <- "In dissertation only (NOT in baseline)"
    extra_rows$baseline_significant <- FALSE

    if ("FDR_disser" %in% colnames(extra_rows)) {
      extra_rows$disser_significant <-
        (abs(extra_rows$logFC_disser) > 1) &
        (extra_rows$FDR_disser < 0.05)
    }

    # Add direction columns
    extra_rows$direction_baseline <- NA
    extra_rows$direction_disser <- ifelse(extra_rows$logFC_disser > 0, "Up", "Down")
    extra_rows$direction_match <- NA

    # Add change columns if they exist
    if ("neglog10FDR_baseline" %in% colnames(comparison_table)) {
      extra_rows$neglog10FDR_baseline <- NA
      extra_rows$neglog10FDR_disser <- -log10(extra_rows$FDR_disser + 1e-300)
      extra_rows$neglog10FDR_disser <- pmin(extra_rows$neglog10FDR_disser, 50)
      extra_rows$significance <- ifelse(extra_rows$FDR_disser < 0.05, "Significant in one", "Not significant")
      extra_rows$logFC_change <- NA
      extra_rows$FDR_change <- NA
      extra_rows$total_change <- NA
    }

    # Bind to comparison table
    comparison_table <- rbind(comparison_table, extra_rows)
  }

  # Add genes significant in baseline but not in dissertation
  if (length(genes_baseline_sig_not_in_disser) > 0) {
    baseline_for_merge_extra <- baseline_difexp[baseline_difexp[[gene_col_baseline]] %in% genes_baseline_sig_not_in_disser, ]

    # Create rows for these genes
    if (has_string_id) {
      extra_rows_baseline <- data.frame(
        gene = baseline_for_merge_extra[[gene_col_baseline]],
        STRING_id = baseline_for_merge_extra$STRING_id,
        logFC_baseline = baseline_for_merge_extra$logFC,
        FDR_baseline = baseline_for_merge_extra$adj.P.Val,
        logFC_disser = NA,
        FDR_disser = NA,
        stringsAsFactors = FALSE
      )
    } else {
      extra_rows_baseline <- data.frame(
        gene = baseline_for_merge_extra[[gene_col_baseline]],
        logFC_baseline = baseline_for_merge_extra$logFC,
        FDR_baseline = baseline_for_merge_extra$adj.P.Val,
        logFC_disser = NA,
        FDR_disser = NA,
        stringsAsFactors = FALSE
      )
    }

    extra_rows_baseline$classification <- "Significant in baseline (|logFC| > 1, adj.P.Val < 0.05) but NOT in dissertation"
    extra_rows_baseline$baseline_significant <- TRUE
    extra_rows_baseline$disser_significant <- FALSE

    # Add direction columns
    extra_rows_baseline$direction_baseline <- ifelse(extra_rows_baseline$logFC_baseline > 0, "Up", "Down")
    extra_rows_baseline$direction_disser <- NA
    extra_rows_baseline$direction_match <- NA

    # Add change columns if they exist
    if ("neglog10FDR_baseline" %in% colnames(comparison_table)) {
      extra_rows_baseline$neglog10FDR_baseline <- -log10(extra_rows_baseline$FDR_baseline + 1e-300)
      extra_rows_baseline$neglog10FDR_baseline <- pmin(extra_rows_baseline$neglog10FDR_baseline, 50)
      extra_rows_baseline$neglog10FDR_disser <- NA
      extra_rows_baseline$significance <- "Significant in one"
      extra_rows_baseline$logFC_change <- NA
      extra_rows_baseline$FDR_change <- NA
      extra_rows_baseline$total_change <- NA
    }

    # Bind to comparison table
    comparison_table <- rbind(comparison_table, extra_rows_baseline)
  }

  # Sort by classification and then by significance
  comparison_table <- comparison_table[order(comparison_table$classification,
                                              -comparison_table$baseline_significant,
                                              comparison_table$FDR_baseline), ]

  # Save comprehensive comparison table
  comprehensive_file <- file.path(output_dir, "comprehensive_comparison_table.csv")
  write.csv(comparison_table, comprehensive_file, row.names = FALSE)
  cat("✓ Saved comprehensive comparison table to:", comprehensive_file, "\n")

  # Print summary statistics
  cat("\nComprehensive comparison statistics:\n")
  cat("  Total genes in comparison table:", nrow(comparison_table), "\n")
  classification_counts <- table(comparison_table$classification)
  for (i in seq_along(classification_counts)) {
    cat("    -", names(classification_counts)[i], ":", classification_counts[i], "\n")
  }
  cat("\n")

  # Create a summary of the classifications
  classification_summary <- data.frame(
    Classification = names(classification_counts),
    Count = as.numeric(classification_counts),
    stringsAsFactors = FALSE
  )

  classification_summary_file <- file.path(output_dir, "classification_summary.csv")
  write.csv(classification_summary, classification_summary_file, row.names = FALSE)
  cat("✓ Saved classification summary to:", classification_summary_file, "\n\n")

  # Save specific gene lists
  if (length(genes_in_disser_not_in_baseline) > 0) {
    disser_not_baseline_file <- file.path(output_dir, "genes_in_disser_NOT_in_baseline.txt")
    writeLines(as.character(genes_in_disser_not_in_baseline), disser_not_baseline_file)
    cat("✓ Saved genes in dissertation but NOT in baseline to:", disser_not_baseline_file, "\n")
  }

  if (length(genes_baseline_sig_not_in_disser) > 0) {
    baseline_sig_not_disser_file <- file.path(output_dir, "genes_baseline_significant_NOT_in_disser.txt")
    writeLines(as.character(genes_baseline_sig_not_in_disser), baseline_sig_not_disser_file)
    cat("✓ Saved genes significant in baseline but NOT in dissertation to:", baseline_sig_not_disser_file, "\n")
  }

  cat("\n")
}

# ===================================================================
# STEP 8: SAVE GENE LISTS
# ===================================================================

cat("Saving gene lists...\n")

# Save common genes
common_genes_file <- file.path(output_dir, "common_genes.txt")
writeLines(as.character(genes_common), common_genes_file)
cat("✓ Saved common genes to:", common_genes_file, "\n")

# Save baseline-only genes
baseline_only_file <- file.path(output_dir, "baseline_only_genes.txt")
writeLines(as.character(genes_baseline_only), baseline_only_file)
cat("✓ Saved baseline-only genes to:", baseline_only_file, "\n")

# Save dissertation-only genes
disser_only_file <- file.path(output_dir, "dissertation_only_genes.txt")
writeLines(as.character(genes_disser_only), disser_only_file)
cat("✓ Saved dissertation-only genes to:", disser_only_file, "\n\n")

# ===================================================================
# STEP 9: SUMMARY REPORT
# ===================================================================

cat("Creating summary report...\n")

summary_report <- data.frame(
  Metric = c(
    "Baseline all genes (unfiltered)",
    "Baseline DEGs (|logFC| > 1, adj.P.Val < 0.05)",
    "Dissertation DEGs",
    "Common genes (baseline all vs disser)",
    "Baseline only",
    "Dissertation only",
    "Overlap (% of baseline all)",
    "Overlap (% of dissertation)",
    if (exists("genes_in_disser_not_in_baseline")) "Genes in disser NOT in baseline" else NULL,
    if (exists("genes_baseline_sig_not_in_disser")) "Genes baseline sig NOT in disser" else NULL,
    if (exists("cor_pearson")) "Pearson correlation" else NULL,
    if (exists("cor_spearman")) "Spearman correlation" else NULL,
    if (exists("pct_match")) "Direction concordance (%)" else NULL
  ),
  Value = c(
    length(baseline_genes),
    length(baseline_genes_filtered),
    length(disser_genes),
    length(genes_common),
    length(genes_baseline_only),
    length(genes_disser_only),
    pct_overlap_baseline,
    pct_overlap_disser,
    if (exists("genes_in_disser_not_in_baseline")) length(genes_in_disser_not_in_baseline) else NULL,
    if (exists("genes_baseline_sig_not_in_disser")) length(genes_baseline_sig_not_in_disser) else NULL,
    if (exists("cor_pearson")) round(cor_pearson, 3) else NULL,
    if (exists("cor_spearman")) round(cor_spearman, 3) else NULL,
    if (exists("pct_match")) pct_match else NULL
  ),
  stringsAsFactors = FALSE
)

# Save summary
summary_file <- file.path(output_dir, "comparison_summary.csv")
write.csv(summary_report, summary_file, row.names = FALSE)

cat("\n")
cat("==================================================================\n")
cat("                    COMPARISON SUMMARY                            \n")
cat("==================================================================\n\n")

print(summary_report, row.names = FALSE)

cat("\n✓ Saved summary to:", summary_file, "\n")

cat("\n")
cat("==================================================================\n")
cat("                    COMPARISON COMPLETE                           \n")
cat("==================================================================\n\n")

cat("All output files saved to:", output_dir, "\n\n")
