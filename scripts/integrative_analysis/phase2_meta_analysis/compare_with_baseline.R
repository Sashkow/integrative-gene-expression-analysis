#!/usr/bin/env Rscript
#' Compare Phase 2 Meta-Analysis Results with ComBat+limma Baseline
#'
#' This script compares the meta-analysis results (Metafor, DExMA, RankProd)
#' with the traditional ComBat+limma approach used as the baseline.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(VennDiagram)
})

# ============================================================================
# Configuration
# ============================================================================

baseline_dir <- "output/trim_2_3"
phase2_dir <- "output/phase2_meta"
output_dir <- "output/phase2_meta/comparison"

# Create output directory
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Thresholds
fdr_threshold <- 0.05
logfc_threshold <- 1.0  # abs(logFC) > 1 for baseline comparison

# ============================================================================
# Load Data
# ============================================================================

cat("Loading data...\n")

# Baseline ComBat+limma results (all genes)
baseline_all <- read.csv(file.path(baseline_dir, "difexp", "difexp_all.csv"),
                         stringsAsFactors = FALSE)
cat("  Baseline all: ", nrow(baseline_all), " genes\n")

# Baseline filtered (FDR < 0.05 and |logFC| > 1)
baseline_filtered <- read.csv(file.path(baseline_dir, "difexp", "difexp_filtered.csv"),
                              stringsAsFactors = FALSE)
cat("  Baseline filtered (FDR<0.05, |logFC|>1): ", nrow(baseline_filtered), " genes\n")

# Phase 2 results
metafor_results <- read.delim(file.path(phase2_dir, "metafor_results.tsv"),
                              stringsAsFactors = FALSE)
cat("  Metafor all: ", nrow(metafor_results), " genes\n")

dexma_results <- read.delim(file.path(phase2_dir, "dexma_results.tsv"),
                            stringsAsFactors = FALSE)
cat("  DExMA all: ", nrow(dexma_results), " genes\n")

rankprod_results <- read.delim(file.path(phase2_dir, "rankprod_results.tsv"),
                               stringsAsFactors = FALSE)
cat("  RankProd all: ", nrow(rankprod_results), " genes\n")

high_confidence <- read.delim(file.path(phase2_dir, "high_confidence_degs.tsv"),
                              stringsAsFactors = FALSE)
cat("  High confidence (both DExMA+RankProd): ", nrow(high_confidence), " genes\n")

# ============================================================================
# Standardize Gene IDs
# ============================================================================

# Baseline uses ENTREZID column, Phase 2 uses gene_id (also ENTREZID)
baseline_all$gene_id <- as.character(baseline_all$ENTREZID)
baseline_filtered$gene_id <- as.character(baseline_filtered$ENTREZID)

# ============================================================================
# Define Significant Gene Sets
# ============================================================================

cat("\n=== Defining Significant Gene Sets ===\n")

# Baseline: FDR < 0.05 only
baseline_sig_fdr <- baseline_all$gene_id[baseline_all$adj.P.Val < fdr_threshold]
cat("Baseline FDR<0.05: ", length(baseline_sig_fdr), " genes\n")

# Baseline: FDR < 0.05 AND |logFC| > 1
baseline_sig_strict <- baseline_filtered$gene_id
cat("Baseline FDR<0.05 + |logFC|>1: ", length(baseline_sig_strict), " genes\n")

# Metafor: FDR < 0.05
metafor_sig <- metafor_results$gene_id[metafor_results$fdr < fdr_threshold]
cat("Metafor FDR<0.05: ", length(metafor_sig), " genes\n")

# DExMA: FDR < 0.05
dexma_sig <- dexma_results$gene_id[dexma_results$fdr < fdr_threshold]
cat("DExMA FDR<0.05: ", length(dexma_sig), " genes\n")

# RankProd: PFP < 0.05
rankprod_sig <- rankprod_results$gene_id[rankprod_results$pfp < fdr_threshold]
cat("RankProd PFP<0.05: ", length(rankprod_sig), " genes\n")

# High confidence
high_conf_sig <- high_confidence$gene_id
cat("High confidence: ", length(high_conf_sig), " genes\n")

# ============================================================================
# Apply Same logFC Filter to Meta-Analysis
# ============================================================================

cat("\n=== With |logFC| > 1 Filter ===\n")

# Metafor with logFC filter
metafor_sig_strict <- metafor_results$gene_id[
  metafor_results$fdr < fdr_threshold & abs(metafor_results$logFC) > logfc_threshold
]
cat("Metafor FDR<0.05 + |logFC|>1: ", length(metafor_sig_strict), " genes\n")

# DExMA with logFC filter
dexma_sig_strict <- dexma_results$gene_id[
  dexma_results$fdr < fdr_threshold & abs(dexma_results$logFC) > logfc_threshold
]
cat("DExMA FDR<0.05 + |logFC|>1: ", length(dexma_sig_strict), " genes\n")

# High confidence with logFC filter
high_conf_sig_strict <- high_confidence$gene_id[
  abs(high_confidence$logFC) > logfc_threshold
]
cat("High confidence |logFC|>1: ", length(high_conf_sig_strict), " genes\n")

# ============================================================================
# Overlap Analysis (FDR only)
# ============================================================================

cat("\n=== Overlap Analysis (FDR < 0.05 only) ===\n")

# Baseline vs Metafor
overlap_baseline_metafor <- length(intersect(baseline_sig_fdr, metafor_sig))
cat("Baseline ∩ Metafor: ", overlap_baseline_metafor,
    " (", round(100 * overlap_baseline_metafor / length(baseline_sig_fdr), 1),
    "% of baseline)\n")

# Baseline vs DExMA
overlap_baseline_dexma <- length(intersect(baseline_sig_fdr, dexma_sig))
cat("Baseline ∩ DExMA: ", overlap_baseline_dexma,
    " (", round(100 * overlap_baseline_dexma / length(baseline_sig_fdr), 1),
    "% of baseline)\n")

# Baseline vs RankProd
overlap_baseline_rankprod <- length(intersect(baseline_sig_fdr, rankprod_sig))
cat("Baseline ∩ RankProd: ", overlap_baseline_rankprod,
    " (", round(100 * overlap_baseline_rankprod / length(baseline_sig_fdr), 1),
    "% of baseline)\n")

# Baseline vs High Confidence
overlap_baseline_highconf <- length(intersect(baseline_sig_fdr, high_conf_sig))
cat("Baseline ∩ High Confidence: ", overlap_baseline_highconf,
    " (", round(100 * overlap_baseline_highconf / length(baseline_sig_fdr), 1),
    "% of baseline)\n")

# All three overlap
overlap_all_three <- length(intersect(intersect(baseline_sig_fdr, metafor_sig), dexma_sig))
cat("Baseline ∩ Metafor ∩ DExMA: ", overlap_all_three, " genes\n")

# ============================================================================
# Overlap Analysis (FDR + logFC)
# ============================================================================

cat("\n=== Overlap Analysis (FDR < 0.05 + |logFC| > 1) ===\n")

# Baseline strict vs Metafor strict
overlap_strict_metafor <- length(intersect(baseline_sig_strict, metafor_sig_strict))
cat("Baseline ∩ Metafor (strict): ", overlap_strict_metafor,
    " (", round(100 * overlap_strict_metafor / length(baseline_sig_strict), 1),
    "% of baseline)\n")

# Baseline strict vs DExMA strict
overlap_strict_dexma <- length(intersect(baseline_sig_strict, dexma_sig_strict))
cat("Baseline ∩ DExMA (strict): ", overlap_strict_dexma,
    " (", round(100 * overlap_strict_dexma / length(baseline_sig_strict), 1),
    "% of baseline)\n")

# Baseline strict vs High Confidence strict
overlap_strict_highconf <- length(intersect(baseline_sig_strict, high_conf_sig_strict))
cat("Baseline ∩ High Confidence (strict): ", overlap_strict_highconf,
    " (", round(100 * overlap_strict_highconf / length(baseline_sig_strict), 1),
    "% of baseline)\n")

# ============================================================================
# Direction Agreement Analysis
# ============================================================================

cat("\n=== Direction Agreement Analysis ===\n")

# For genes in common between baseline and DExMA, check direction agreement
common_baseline_dexma <- intersect(baseline_sig_fdr, dexma_sig)

if (length(common_baseline_dexma) > 0) {
  baseline_subset <- baseline_all[baseline_all$gene_id %in% common_baseline_dexma, ]
  dexma_subset <- dexma_results[dexma_results$gene_id %in% common_baseline_dexma, ]

  # Merge on gene_id
  merged <- merge(baseline_subset[, c("gene_id", "logFC")],
                  dexma_subset[, c("gene_id", "logFC")],
                  by = "gene_id", suffixes = c("_baseline", "_dexma"))

  # Check direction agreement (same sign)
  merged$same_direction <- sign(merged$logFC_baseline) == sign(merged$logFC_dexma)
  direction_agree <- sum(merged$same_direction)

  cat("Baseline vs DExMA direction agreement: ", direction_agree, "/",
      nrow(merged), " (", round(100 * direction_agree / nrow(merged), 1), "%)\n")
}

# For genes in common between baseline and Metafor
common_baseline_metafor <- intersect(baseline_sig_fdr, metafor_sig)

if (length(common_baseline_metafor) > 0) {
  baseline_subset <- baseline_all[baseline_all$gene_id %in% common_baseline_metafor, ]
  metafor_subset <- metafor_results[metafor_results$gene_id %in% common_baseline_metafor, ]

  merged <- merge(baseline_subset[, c("gene_id", "logFC")],
                  metafor_subset[, c("gene_id", "logFC")],
                  by = "gene_id", suffixes = c("_baseline", "_metafor"))

  merged$same_direction <- sign(merged$logFC_baseline) == sign(merged$logFC_metafor)
  direction_agree <- sum(merged$same_direction)

  cat("Baseline vs Metafor direction agreement: ", direction_agree, "/",
      nrow(merged), " (", round(100 * direction_agree / nrow(merged), 1), "%)\n")
}

# ============================================================================
# Correlation of Effect Sizes
# ============================================================================

cat("\n=== Effect Size Correlation ===\n")

# Common genes between baseline and DExMA (all genes, not just significant)
common_all <- intersect(baseline_all$gene_id, dexma_results$gene_id)
cat("Genes in both baseline and DExMA: ", length(common_all), "\n")

baseline_common <- baseline_all[match(common_all, baseline_all$gene_id), ]
dexma_common <- dexma_results[match(common_all, dexma_results$gene_id), ]

# Correlation of logFC
cor_logfc <- cor(baseline_common$logFC, dexma_common$logFC, use = "complete.obs")
cat("Correlation of logFC (baseline vs DExMA): ", round(cor_logfc, 3), "\n")

# Correlation with Metafor
common_metafor <- intersect(baseline_all$gene_id, metafor_results$gene_id)
baseline_mf <- baseline_all[match(common_metafor, baseline_all$gene_id), ]
metafor_mf <- metafor_results[match(common_metafor, metafor_results$gene_id), ]

cor_logfc_mf <- cor(baseline_mf$logFC, metafor_mf$logFC, use = "complete.obs")
cat("Correlation of logFC (baseline vs Metafor): ", round(cor_logfc_mf, 3), "\n")

# ============================================================================
# Unique Discoveries by Each Method
# ============================================================================

cat("\n=== Unique Discoveries ===\n")

# Genes found only by meta-analysis (not in baseline FDR < 0.05)
unique_dexma <- setdiff(dexma_sig, baseline_sig_fdr)
unique_metafor <- setdiff(metafor_sig, baseline_sig_fdr)
unique_highconf <- setdiff(high_conf_sig, baseline_sig_fdr)

cat("Unique to DExMA (not in baseline): ", length(unique_dexma), " genes\n")
cat("Unique to Metafor (not in baseline): ", length(unique_metafor), " genes\n")
cat("Unique to High Confidence (not in baseline): ", length(unique_highconf), " genes\n")

# Genes found only by baseline (not in meta-analysis)
unique_baseline_vs_dexma <- setdiff(baseline_sig_fdr, dexma_sig)
unique_baseline_vs_metafor <- setdiff(baseline_sig_fdr, metafor_sig)
unique_baseline_vs_highconf <- setdiff(baseline_sig_fdr, high_conf_sig)

cat("Unique to baseline (not in DExMA): ", length(unique_baseline_vs_dexma), " genes\n")
cat("Unique to baseline (not in Metafor): ", length(unique_baseline_vs_metafor), " genes\n")
cat("Unique to baseline (not in High Conf): ", length(unique_baseline_vs_highconf), " genes\n")

# ============================================================================
# Generate Plots
# ============================================================================

cat("\n=== Generating Plots ===\n")

# 1. Scatter plot of logFC values
pdf(file.path(output_dir, "logfc_correlation_baseline_dexma.pdf"), width = 8, height = 8)
plot(baseline_common$logFC, dexma_common$logFC,
     xlab = "logFC (ComBat+limma baseline)",
     ylab = "logFC (DExMA meta-analysis)",
     main = paste0("Effect Size Correlation\nr = ", round(cor_logfc, 3)),
     pch = 16, col = rgb(0, 0, 0, 0.3), cex = 0.5)
abline(0, 1, col = "red", lty = 2)
abline(h = 0, v = 0, col = "gray", lty = 3)
dev.off()
cat("  Saved logfc_correlation_baseline_dexma.pdf\n")

# 2. Scatter plot baseline vs Metafor
pdf(file.path(output_dir, "logfc_correlation_baseline_metafor.pdf"), width = 8, height = 8)
plot(baseline_mf$logFC, metafor_mf$logFC,
     xlab = "logFC (ComBat+limma baseline)",
     ylab = "logFC (Metafor meta-analysis)",
     main = paste0("Effect Size Correlation\nr = ", round(cor_logfc_mf, 3)),
     pch = 16, col = rgb(0, 0, 0, 0.3), cex = 0.5)
abline(0, 1, col = "red", lty = 2)
abline(h = 0, v = 0, col = "gray", lty = 3)
dev.off()
cat("  Saved logfc_correlation_baseline_metafor.pdf\n")

# 3. Venn diagram (FDR only)
if (requireNamespace("VennDiagram", quietly = TRUE)) {
  # Suppress log file creation
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

  venn_data <- list(
    "ComBat+limma" = baseline_sig_fdr,
    "DExMA" = dexma_sig,
    "Metafor" = metafor_sig
  )

  pdf(file.path(output_dir, "venn_three_methods.pdf"), width = 10, height = 10)
  grid.newpage()
  venn_plot <- venn.diagram(
    x = venn_data,
    filename = NULL,
    fill = c("lightblue", "lightgreen", "lightyellow"),
    alpha = 0.5,
    cat.cex = 1.2,
    cex = 1.5,
    main = "Overlap of Significant DEGs (FDR < 0.05)",
    main.cex = 1.5
  )
  grid.draw(venn_plot)
  dev.off()
  cat("  Saved venn_three_methods.pdf\n")
}

# ============================================================================
# Summary Statistics Table
# ============================================================================

cat("\n=== Summary Statistics ===\n")

summary_df <- data.frame(
  Method = c("ComBat+limma (baseline)", "Metafor", "DExMA", "RankProd", "High Confidence"),
  Total_Genes = c(nrow(baseline_all), nrow(metafor_results), nrow(dexma_results),
                  nrow(rankprod_results), nrow(high_confidence)),
  Significant_FDR = c(length(baseline_sig_fdr), length(metafor_sig), length(dexma_sig),
                      length(rankprod_sig), length(high_conf_sig)),
  Significant_Strict = c(length(baseline_sig_strict), length(metafor_sig_strict),
                         length(dexma_sig_strict), NA, length(high_conf_sig_strict)),
  Overlap_with_Baseline_FDR = c(NA, overlap_baseline_metafor, overlap_baseline_dexma,
                                overlap_baseline_rankprod, overlap_baseline_highconf),
  Pct_Baseline_Recovered = c(NA,
                              round(100 * overlap_baseline_metafor / length(baseline_sig_fdr), 1),
                              round(100 * overlap_baseline_dexma / length(baseline_sig_fdr), 1),
                              round(100 * overlap_baseline_rankprod / length(baseline_sig_fdr), 1),
                              round(100 * overlap_baseline_highconf / length(baseline_sig_fdr), 1))
)

print(summary_df)

# Save summary table
write.csv(summary_df, file.path(output_dir, "comparison_summary.csv"), row.names = FALSE)
cat("\nSaved comparison_summary.csv\n")

# ============================================================================
# Detailed Comparison Report
# ============================================================================

report_lines <- c(
  "Phase 2 Meta-Analysis vs ComBat+limma Baseline Comparison",
  "=========================================================",
  paste0("Generated: ", Sys.time()),
  "",
  "OVERVIEW",
  "--------",
  paste0("Baseline (ComBat+limma) genes: ", nrow(baseline_all)),
  paste0("Baseline significant (FDR<0.05): ", length(baseline_sig_fdr)),
  paste0("Baseline strict (FDR<0.05 + |logFC|>1): ", length(baseline_sig_strict)),
  "",
  "META-ANALYSIS RESULTS",
  "--------------------",
  paste0("Metafor significant (FDR<0.05): ", length(metafor_sig)),
  paste0("DExMA significant (FDR<0.05): ", length(dexma_sig)),
  paste0("RankProd significant (PFP<0.05): ", length(rankprod_sig)),
  paste0("High Confidence (DExMA + RankProd): ", length(high_conf_sig)),
  "",
  "OVERLAP WITH BASELINE",
  "--------------------",
  paste0("Baseline genes recovered by Metafor: ", overlap_baseline_metafor,
         " (", round(100 * overlap_baseline_metafor / length(baseline_sig_fdr), 1), "%)"),
  paste0("Baseline genes recovered by DExMA: ", overlap_baseline_dexma,
         " (", round(100 * overlap_baseline_dexma / length(baseline_sig_fdr), 1), "%)"),
  paste0("Baseline genes recovered by RankProd: ", overlap_baseline_rankprod,
         " (", round(100 * overlap_baseline_rankprod / length(baseline_sig_fdr), 1), "%)"),
  paste0("Baseline genes recovered by High Confidence: ", overlap_baseline_highconf,
         " (", round(100 * overlap_baseline_highconf / length(baseline_sig_fdr), 1), "%)"),
  "",
  "EFFECT SIZE CORRELATION",
  "-----------------------",
  paste0("Correlation (baseline vs DExMA): ", round(cor_logfc, 3)),
  paste0("Correlation (baseline vs Metafor): ", round(cor_logfc_mf, 3)),
  "",
  "INTERPRETATION",
  "--------------",
  "The meta-analysis methods show different characteristics compared to the",
  "ComBat+limma baseline:",
  "",
  "- Metafor (effect-size random effects) provides the broadest gene coverage",
  "  by utilizing data from all 7 studies, even when gene measurements vary",
  "  across studies.",
  "",
  "- DExMA and RankProd require balanced studies (samples in both conditions)",
  "  so they use only 2 of 7 studies, resulting in fewer genes but potentially",
  "  more reliable cross-study validation.",
  "",
  "- The high confidence DEGs (significant in both DExMA and RankProd with",
  "  direction agreement) represent the most robust findings.",
  "",
  "- The baseline ComBat+limma approach merges data before analysis, which",
  "  can mask study-specific effects but provides a unified view.",
  ""
)

writeLines(report_lines, file.path(output_dir, "comparison_report.txt"))
cat("\nSaved comparison_report.txt\n")

cat("\n=== Comparison Complete ===\n")
