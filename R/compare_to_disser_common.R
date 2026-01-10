#' Common Comparison Functions for Dissertation Data
#'
#' Shared functions for comparing differential expression results
#' with dissertation data (any sheet).
#'
#' @author Expression Integration Pipeline
#' @date 2025-01-07

library(readxl)
library(ggplot2)
library(VennDiagram)

#' Run full comparison between new results and dissertation data
#'
#' @param results_dir Directory containing difexp results (difexp_final.csv, difexp/difexp_all.csv)
#' @param disser_file Path to dissertation xlsx file
#' @param disser_sheet Sheet number to read from dissertation file (1 for 1_2, 2 for 2_3)
#' @param output_dir Output directory for comparison results
#' @param comparison_name Name for the comparison (used in titles/labels)
run_dissertation_comparison <- function(results_dir,
                                         disser_file,
                                         disser_sheet,
                                         output_dir,
                                         comparison_name) {

  cat("\n=== Comparing", comparison_name, "to Dissertation Data (Sheet", disser_sheet, ") ===\n\n")

  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n\n")

  # ===================================================================
  # STEP 1: LOAD NEW RESULTS
  # ===================================================================

  cat("Loading new results...\n")

  difexp_file <- file.path(results_dir, "difexp_final.csv")
  difexp_all_file <- file.path(results_dir, "difexp", "difexp_all.csv")

  # Fallback for baseline structure

  if (!file.exists(difexp_file)) {
    difexp_file <- file.path(results_dir, "difexp_genes.csv")
  }
  if (!file.exists(difexp_all_file)) {
    difexp_all_file <- file.path(results_dir, "difexp_genes_unfiltered.csv")
  }

  if (!file.exists(difexp_file)) {
    stop("Differential expression file not found in: ", results_dir)
  }

  difexp <- read.csv(difexp_file, stringsAsFactors = FALSE)
  difexp_all <- read.csv(difexp_all_file, stringsAsFactors = FALSE)

  cat("  Loaded new results:\n")
  cat("    Filtered DEGs:", nrow(difexp), "\n")
  cat("    All genes:", nrow(difexp_all), "\n\n")

  # ===================================================================
  # STEP 2: LOAD DISSERTATION DATA
  # ===================================================================

  cat("Loading dissertation data (sheet", disser_sheet, ")...\n")

  if (!file.exists(disser_file)) {
    stop("Dissertation file not found: ", disser_file)
  }

  disser_data <- read_excel(disser_file, sheet = disser_sheet)

  cat("  Loaded dissertation data:\n")
  cat("    Total genes:", nrow(disser_data), "\n")
  cat("    Columns:", paste(colnames(disser_data), collapse = ", "), "\n\n")

  # ===================================================================
  # STEP 3: PREPARE DATA FOR COMPARISON
  # ===================================================================

  cat("Preparing data for comparison...\n")

  # Identify gene column
  gene_col_new <- if ("ENTREZID" %in% colnames(difexp)) "ENTREZID" else "SYMBOL"
  gene_col_disser <- if ("ENTREZID" %in% colnames(disser_data)) "ENTREZID" else
                     if ("SYMBOL" %in% colnames(disser_data)) "SYMBOL" else colnames(disser_data)[1]

  cat("  Using gene column - New:", gene_col_new, ", Dissertation:", gene_col_disser, "\n")

  # Extract gene lists
  new_genes <- difexp_all[[gene_col_new]]
  new_genes_filtered <- difexp[[gene_col_new]]
  disser_genes <- disser_data[[gene_col_disser]]

  # Remove NAs
  new_genes <- new_genes[!is.na(new_genes)]
  new_genes_filtered <- new_genes_filtered[!is.na(new_genes_filtered)]
  disser_genes <- disser_genes[!is.na(disser_genes)]

  cat("  Gene counts: New all =", length(new_genes),
      ", New filtered =", length(new_genes_filtered),
      ", Dissertation =", length(disser_genes), "\n\n")

  # ===================================================================
  # STEP 4: OVERLAP ANALYSIS (DEGs vs DEGs)
  # ===================================================================

  cat("Computing overlaps (filtered DEGs vs dissertation DEGs)...\n")

  # Use filtered DEGs for Venn comparison
  genes_common <- intersect(new_genes_filtered, disser_genes)
  genes_new_only <- setdiff(new_genes_filtered, disser_genes)
  genes_disser_only <- setdiff(disser_genes, new_genes_filtered)

  pct_overlap_new <- round(100 * length(genes_common) / length(new_genes_filtered), 1)
  pct_overlap_disser <- round(100 * length(genes_common) / length(disser_genes), 1)

  cat("  New filtered DEGs:", length(new_genes_filtered), "\n")
  cat("  Dissertation DEGs:", length(disser_genes), "\n")
  cat("  Common:", length(genes_common), "\n")
  cat("  New only:", length(genes_new_only), "\n")
  cat("  Dissertation only:", length(genes_disser_only), "\n")
  cat("  Overlap: ", pct_overlap_disser, "% of dissertation DEGs recovered\n\n")

  # Venn diagram (DEGs vs DEGs)
  venn_plot <- venn.diagram(
    x = list(`New DEGs` = new_genes_filtered, `Disser DEGs` = disser_genes),
    filename = NULL,
    fill = c("#440154ff", "#21908dff"),
    alpha = 0.5, cex = 1.5, cat.cex = 1.5, cat.fontface = "bold",
    cat.default.pos = "outer", cat.dist = c(0.055, 0.055)
  )

  venn_file <- file.path(output_dir, "venn_diagram.png")
  png(venn_file, width = 8, height = 8, units = "in", res = 300)
  grid::grid.draw(venn_plot)
  dev.off()
  cat("  Saved Venn diagram\n")

  # ===================================================================
  # STEP 5: LOGFC COMPARISON
  # ===================================================================

  merged_data <- NULL
  cor_pearson <- NA
  cor_spearman <- NA
  pct_match <- NA

  if ("logFC" %in% colnames(difexp_all) && "logFC" %in% colnames(disser_data)) {

    cat("Comparing logFC values...\n")

    use_string_id <- "STRING_id" %in% colnames(difexp_all) &&
                     "STRING_id" %in% colnames(disser_data)

    if (use_string_id) {
      new_for_merge <- difexp_all[, c(gene_col_new, "STRING_id", "logFC", "adj.P.Val")]
      colnames(new_for_merge) <- c("gene", "STRING_id", "logFC_new", "FDR_new")

      disser_for_merge <- disser_data[, c(gene_col_disser, "STRING_id", "logFC")]
      colnames(disser_for_merge) <- c("gene", "STRING_id", "logFC_disser")

      # Get FDR from dissertation
      if ("adj.P.Val" %in% colnames(disser_data)) {
        disser_for_merge$FDR_disser <- disser_data$adj.P.Val
      }

      merged_data <- merge(new_for_merge, disser_for_merge, by = c("gene", "STRING_id"), all = FALSE)
    } else {
      new_for_merge <- difexp_all[, c(gene_col_new, "logFC", "adj.P.Val")]
      colnames(new_for_merge) <- c("gene", "logFC_new", "FDR_new")
      new_for_merge <- new_for_merge[!duplicated(new_for_merge$gene), ]

      disser_for_merge <- data.frame(
        gene = disser_data[[gene_col_disser]],
        logFC_disser = disser_data$logFC,
        stringsAsFactors = FALSE
      )
      if ("adj.P.Val" %in% colnames(disser_data)) {
        disser_for_merge$FDR_disser <- disser_data$adj.P.Val
      }
      disser_for_merge <- disser_for_merge[!duplicated(disser_for_merge$gene), ]

      merged_data <- merge(new_for_merge, disser_for_merge, by = "gene", all = FALSE)
    }

    cat("  Merged genes:", nrow(merged_data), "\n")

    # Correlation
    cor_pearson <- cor(merged_data$logFC_new, merged_data$logFC_disser, use = "complete.obs")
    cor_spearman <- cor(merged_data$logFC_new, merged_data$logFC_disser,
                        method = "spearman", use = "complete.obs")

    cat("  Pearson r:", round(cor_pearson, 3), "\n")
    cat("  Spearman r:", round(cor_spearman, 3), "\n")

    # Direction concordance
    merged_data$direction_new <- ifelse(merged_data$logFC_new > 0, "Up", "Down")
    merged_data$direction_disser <- ifelse(merged_data$logFC_disser > 0, "Up", "Down")
    merged_data$direction_match <- merged_data$direction_new == merged_data$direction_disser

    n_match <- sum(merged_data$direction_match, na.rm = TRUE)
    pct_match <- round(100 * n_match / nrow(merged_data), 1)
    cat("  Direction concordance:", pct_match, "%\n\n")

    # Scatter plot
    p <- ggplot(merged_data, aes(x = logFC_disser, y = logFC_new)) +
      geom_point(alpha = 0.5, size = 2) +
      geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
      geom_smooth(method = "lm", color = "blue", se = TRUE) +
      labs(
        title = paste("logFC Comparison:", comparison_name, "vs Dissertation"),
        subtitle = paste0("N = ", nrow(merged_data),
                         " | Pearson r = ", round(cor_pearson, 3),
                         " | Spearman r = ", round(cor_spearman, 3)),
        x = "logFC (Dissertation)", y = paste("logFC (", comparison_name, ")")
      ) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
            plot.subtitle = element_text(hjust = 0.5, size = 10))

    scatter_file <- file.path(output_dir, "logFC_scatter.png")
    png(scatter_file, width = 8, height = 8, units = "in", res = 300)
    print(p)
    dev.off()
    cat("  Saved scatter plot\n")

    # Save merged data
    write.csv(merged_data, file.path(output_dir, "common_genes_comparison.csv"), row.names = FALSE)
  }

  # ===================================================================
  # STEP 6: SAVE GENE LISTS
  # ===================================================================

  cat("Saving gene lists...\n")

  writeLines(as.character(genes_common), file.path(output_dir, "common_genes.txt"))
  writeLines(as.character(genes_new_only), file.path(output_dir, "new_only_genes.txt"))
  writeLines(as.character(genes_disser_only), file.path(output_dir, "dissertation_only_genes.txt"))

  genes_new_sig_not_in_disser <- setdiff(new_genes_filtered, disser_genes)
  if (length(genes_new_sig_not_in_disser) > 0) {
    writeLines(as.character(genes_new_sig_not_in_disser),
               file.path(output_dir, "new_significant_NOT_in_disser.txt"))
  }

  genes_in_disser_not_in_new <- setdiff(disser_genes, new_genes_filtered)
  if (length(genes_in_disser_not_in_new) > 0) {
    writeLines(as.character(genes_in_disser_not_in_new),
               file.path(output_dir, "genes_in_disser_NOT_in_new.txt"))
  }

  # ===================================================================
  # STEP 7: SUMMARY REPORT
  # ===================================================================

  summary_report <- data.frame(
    Metric = c(
      "New results - all genes",
      "New results - filtered DEGs",
      "Dissertation DEGs",
      "Common genes",
      "New only",
      "Dissertation only",
      "Overlap (% of new)",
      "Overlap (% of dissertation)",
      "Genes in disser NOT in new",
      "Genes new sig NOT in disser",
      "Pearson correlation",
      "Spearman correlation",
      "Direction concordance (%)"
    ),
    Value = c(
      length(new_genes),
      length(new_genes_filtered),
      length(disser_genes),
      length(genes_common),
      length(genes_new_only),
      length(genes_disser_only),
      pct_overlap_new,
      pct_overlap_disser,
      length(genes_in_disser_not_in_new),
      length(genes_new_sig_not_in_disser),
      round(cor_pearson, 3),
      round(cor_spearman, 3),
      pct_match
    ),
    stringsAsFactors = FALSE
  )

  summary_file <- file.path(output_dir, "comparison_summary.csv")
  write.csv(summary_report, summary_file, row.names = FALSE)

  cat("\n")
  cat("==================================================================\n")
  cat("  ", comparison_name, " vs DISSERTATION COMPARISON SUMMARY\n")
  cat("==================================================================\n\n")
  print(summary_report, row.names = FALSE)
  cat("\nAll output saved to:", output_dir, "\n\n")

  return(summary_report)
}
