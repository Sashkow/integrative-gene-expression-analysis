#!/usr/bin/env Rscript
#
# Phase 5 validation plots (10 plots).
# Reads TSVs from test1, test1b, test2, test3.
#
# Usage: Rscript plot_validation_results.R [--config=config_validation.yaml]

script_dir <- if (length(grep("--file=", commandArgs(FALSE), value = TRUE)) > 0) {
  dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
} else {
  "scripts/integrative_analysis/phase5_validation"
}
source(file.path(script_dir, "subsampling_helpers.R"))

default_config <- file.path(script_dir, "config_validation.yaml")
config <- parse_config_arg(default_config)

input_dir <- config$paths$output
plot_dir <- file.path(input_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

ref_de <- read.delim(config$paths$reference_de)
ref_n_deg <- sum(
  ref_de$adj.P.Val < 0.05 & abs(ref_de$logFC) > 1, na.rm = TRUE
)
ref_n_deg_fdr <- sum(ref_de$adj.P.Val < 0.05, na.rm = TRUE)

bal_n_deg <- NA
bal_de_path <- config$paths$balanced_reference_de
if (!is.null(bal_de_path) && file.exists(bal_de_path)) {
  bal_de <- read.delim(bal_de_path)
  bal_n_deg <- sum(
    bal_de$adj.P.Val < 0.05 & abs(bal_de$logFC) > 1, na.rm = TRUE
  )
  bal_n_deg_fdr <- sum(bal_de$adj.P.Val < 0.05, na.rm = TRUE)
}

# ---- Load all data ----
t1_file <- file.path(input_dir, "test1_first_trim_subsample.tsv")
t1b_file <- file.path(input_dir, "test1b_vs_balanced.tsv")
t2_file <- file.path(input_dir, "test2_balanced_subsample.tsv")
t3_file <- file.path(input_dir, "test3_split_half.tsv")

has_t1 <- file.exists(t1_file)
has_t1b <- file.exists(t1b_file)
has_t2 <- file.exists(t2_file)
has_t3 <- file.exists(t3_file)

if (has_t1) {
  t1 <- read.delim(t1_file, stringsAsFactors = FALSE)
  t1_ok <- t1[t1$status == "ok", ]
}
if (has_t1b) {
  t1b <- read.delim(t1b_file, stringsAsFactors = FALSE)
  t1b_ok <- t1b[t1b$status == "ok", ]
}
if (has_t2) {
  t2 <- read.delim(t2_file, stringsAsFactors = FALSE)
  t2_ok <- t2[t2$status == "ok", ]
}
if (has_t3) {
  t3 <- read.delim(t3_file, stringsAsFactors = FALSE)
  t3_ok <- t3[t3$status_a == "ok" & t3$status_b == "ok", ]
}

# ============================================================
# 01: Test 1 — DEG count by first-trimester sample size
# ============================================================
if (has_t1) {
  png(file.path(plot_dir, "01_test1_deg_count.png"),
      width = 1200, height = 800, res = 150)
  boxplot(n_deg ~ N_1st, data = t1_ok, col = "lightblue",
          xlab = "N first-trimester samples",
          ylab = "Number of DEGs",
          main = "DEG count vs first-trimester sample size")
  abline(h = ref_n_deg, lty = 2, col = "red", lwd = 2)
  abline(h = bal_n_deg, lty = 3, col = "blue", lwd = 2)
  legend("topright",
         legend = c(paste0("Full run (", ref_n_deg, ")"),
                    paste0("Balanced (", bal_n_deg, ")")),
         lty = c(2, 3), col = c("red", "blue"), lwd = 2, bty = "n")
  dev.off()
  cat("01 saved\n")
}

# ============================================================
# 02: Test 1 — Jaccard vs full by N_1st
# ============================================================
if (has_t1) {
  png(file.path(plot_dir, "02_test1_jaccard_vs_full.png"),
      width = 1200, height = 800, res = 150)
  boxplot(jaccard_vs_full ~ N_1st, data = t1_ok, col = "lightgreen",
          xlab = "N first-trimester samples",
          ylab = "Jaccard index vs full run",
          main = "DEG overlap with full result",
          ylim = c(0, 1))
  dev.off()
  cat("02 saved\n")
}

# ============================================================
# 03: Retention — 4-panel
#   full DEGs, balanced DEGs, full FDR-only, balanced FDR-only
# ============================================================
if (has_t1 && has_t1b) {
  png(file.path(plot_dir, "03_retention_all.png"),
      width = 2400, height = 800, res = 150)
  par(mfrow = c(1, 4), mar = c(5, 4, 3, 1))

  boxplot(overlap_vs_full ~ N_1st, data = t1_ok,
          col = "lightblue",
          xlab = "N first-trimester samples",
          ylab = "Retention",
          main = paste0("Full DEGs (", ref_n_deg, ")"),
          ylim = c(0, 1))

  boxplot(overlap_vs_balanced ~ N_1st, data = t1b_ok,
          col = "lightsalmon",
          xlab = "N first-trimester samples",
          ylab = "Retention",
          main = paste0("Balanced DEGs (", bal_n_deg, ")"),
          ylim = c(0, 1))

  boxplot(overlap_fdr_only ~ N_1st, data = t1_ok,
          col = "lightblue",
          xlab = "N first-trimester samples",
          ylab = "Retention",
          main = paste0("Full FDR-only (", ref_n_deg_fdr, ")"),
          ylim = c(0, 1))

  if ("overlap_fdr_vs_balanced" %in% colnames(t1b_ok)) {
    boxplot(overlap_fdr_vs_balanced ~ N_1st, data = t1b_ok,
            col = "lightsalmon",
            xlab = "N first-trimester samples",
            ylab = "Retention",
            main = paste0("Balanced FDR-only (",
                           bal_n_deg_fdr, ")"),
            ylim = c(0, 1))
  } else {
    plot.new()
    text(0.5, 0.5, "overlap_fdr_vs_balanced\nnot available",
         cex = 1.2)
  }

  dev.off()
  cat("03 saved\n")
}

# ============================================================
# 04: logFC CCC — 2-panel (vs full + vs balanced)
# ============================================================
if (has_t1b) {
  ccc_full_col <- if ("logfc_ccc_vs_full" %in% colnames(t1b_ok))
    "logfc_ccc_vs_full" else "logfc_r_vs_full"
  ccc_bal_col <- if ("logfc_ccc_vs_balanced" %in% colnames(t1b_ok))
    "logfc_ccc_vs_balanced" else "logfc_r_vs_balanced"
  ccc_label <- if (grepl("ccc", ccc_full_col)) "Lin's CCC" else
    "Pearson r"

  png(file.path(plot_dir, "04_logfc_ccc.png"),
      width = 1800, height = 800, res = 150)
  par(mfrow = c(1, 2), mar = c(5, 4, 3, 1))

  boxplot(t1b_ok[[ccc_full_col]] ~ t1b_ok$N_1st,
          col = "lightblue",
          xlab = "N first-trimester samples",
          ylab = paste(ccc_label, "(logFC)"),
          main = "logFC agreement vs full",
          ylim = c(0, 1))

  boxplot(t1b_ok[[ccc_bal_col]] ~ t1b_ok$N_1st,
          col = "lightsalmon",
          xlab = "N first-trimester samples",
          ylab = paste(ccc_label, "(logFC)"),
          main = "logFC agreement vs balanced",
          ylim = c(0, 1))

  dev.off()
  cat("04 saved\n")
}

# ============================================================
# 05: FDR+logFC vs FDR-only Jaccard convergence
# ============================================================
if (has_t1 && "jaccard_fdr_only" %in% colnames(t1_ok)) {
  t1_med <- aggregate(jaccard_vs_full ~ N_1st,
                       data = t1_ok, FUN = median)
  t1_med_fdr <- aggregate(jaccard_fdr_only ~ N_1st,
                           data = t1_ok, FUN = median)

  png(file.path(plot_dir, "05_fdr_vs_logfc_stability.png"),
      width = 1200, height = 800, res = 150)
  plot(t1_med$N_1st, t1_med$jaccard_vs_full,
       type = "b", pch = 19, col = "blue", ylim = c(0, 1),
       xlab = "N first-trimester samples",
       ylab = "Median Jaccard vs full",
       main = "DEG list stability: FDR+logFC vs FDR-only")
  points(t1_med_fdr$N_1st, t1_med_fdr$jaccard_fdr_only,
         type = "b", pch = 17, col = "red")
  legend("bottomright",
         legend = c(paste0("FDR + |logFC|>1 (", ref_n_deg, " ref)"),
                    paste0("FDR only (", ref_n_deg_fdr, " ref)")),
         col = c("blue", "red"), pch = c(19, 17), bty = "n")
  dev.off()
  cat("05 saved\n")
}

# ============================================================
# 06: Test 2 — balanced subsampling (3-panel)
# ============================================================
if (has_t2) {
  png(file.path(plot_dir, "06_test2_balanced_subsampling.png"),
      width = 1800, height = 600, res = 150)
  par(mfrow = c(1, 3), mar = c(5, 4, 3, 1))

  boxplot(n_deg ~ n_per_trim, data = t2_ok, col = "lightsalmon",
          xlab = "N per trimester", ylab = "Number of DEGs",
          main = "Balanced: DEG count")
  abline(h = ref_n_deg, lty = 2, col = "red", lwd = 2)

  boxplot(jaccard_vs_full ~ n_per_trim, data = t2_ok,
          col = "lightsalmon",
          xlab = "N per trimester", ylab = "Jaccard vs full",
          main = "Balanced: overlap", ylim = c(0, 1))

  ccc_col <- if ("logfc_ccc_vs_full" %in% colnames(t2_ok))
    "logfc_ccc_vs_full" else "logfc_pearson_vs_full"
  ccc_lab <- if (grepl("ccc", ccc_col)) "Lin's CCC" else
    "Pearson r"
  boxplot(t2_ok[[ccc_col]] ~ t2_ok$n_per_trim,
          col = "lightsalmon",
          xlab = "N per trimester",
          ylab = paste(ccc_lab, "(logFC)"),
          main = "Balanced: logFC agreement", ylim = c(0, 1))

  dev.off()
  cat("06 saved\n")
}

# ============================================================
# 07: Test 3 — split-half (2x2)
# ============================================================
if (has_t3) {
  png(file.path(plot_dir, "07_test3_split_half.png"),
      width = 1600, height = 1200, res = 150)
  par(mfrow = c(2, 2), mar = c(5, 4, 3, 1))

  # Top-left: Jaccard A vs B
  hist(t3_ok$jaccard_a_vs_b, breaks = 20, col = "plum",
       xlab = "Jaccard (half A vs half B)",
       main = "DEG overlap between halves", xlim = c(0, 1))
  abline(v = median(t3_ok$jaccard_a_vs_b, na.rm = TRUE),
         lty = 2, col = "red", lwd = 2)

  # Top-right: CCC A vs B
  if ("logfc_ccc_a_vs_b" %in% colnames(t3_ok)) {
    hist(t3_ok$logfc_ccc_a_vs_b, breaks = 20, col = "plum",
         xlab = "Lin's CCC (half A vs half B)",
         main = "logFC agreement between halves", xlim = c(0, 1))
    abline(v = median(t3_ok$logfc_ccc_a_vs_b, na.rm = TRUE),
           lty = 2, col = "red", lwd = 2)
  } else {
    hist(t3_ok$logfc_pearson_a_vs_b, breaks = 20, col = "plum",
         xlab = "Pearson r (half A vs half B)",
         main = "logFC correlation between halves", xlim = c(0, 1))
    abline(v = median(t3_ok$logfc_pearson_a_vs_b, na.rm = TRUE),
           lty = 2, col = "red", lwd = 2)
  }

  # Bottom-left: retention of full DEGs per half
  if ("overlap_a_vs_full" %in% colnames(t3_ok)) {
    boxplot(
      list("Half A" = t3_ok$overlap_a_vs_full,
           "Half B" = t3_ok$overlap_b_vs_full,
           "Both" = c(t3_ok$overlap_a_vs_full,
                      t3_ok$overlap_b_vs_full)),
      col = c("plum", "plum", "lightyellow"),
      ylab = "Retention of full-run DEGs",
      main = "FDR+logFC retention per half", ylim = c(0, 1))
  } else {
    all_j <- c(t3_ok$jaccard_a_vs_full, t3_ok$jaccard_b_vs_full)
    hist(all_j, breaks = 20, col = "plum",
         xlab = "Jaccard (half vs full)",
         main = "Each half vs full run", xlim = c(0, 1))
    abline(v = median(all_j, na.rm = TRUE),
           lty = 2, col = "red", lwd = 2)
  }

  # Bottom-right: FDR-only retention per half
  if ("overlap_fdr_a_vs_full" %in% colnames(t3_ok)) {
    boxplot(
      list("Half A" = t3_ok$overlap_fdr_a_vs_full,
           "Half B" = t3_ok$overlap_fdr_b_vs_full,
           "Both" = c(t3_ok$overlap_fdr_a_vs_full,
                      t3_ok$overlap_fdr_b_vs_full)),
      col = c("plum", "plum", "lightyellow"),
      ylab = "Retention of full-run FDR-only DEGs",
      main = "FDR-only retention per half", ylim = c(0, 1))
  } else {
    boxplot(
      list("FDR+logFC\nA vs B" = t3_ok$jaccard_a_vs_b,
           "FDR only\nA vs B" = t3_ok$jaccard_fdr_a_vs_b),
      col = c("lightyellow", "plum"),
      main = "FDR+logFC vs FDR-only", ylim = c(0, 1),
      ylab = "Jaccard (A vs B)")
  }

  dev.off()
  cat("07 saved\n")
}

# ============================================================
# 08: Test 3 — FDR+logFC vs FDR-only comparison
# ============================================================
if (has_t3 && "jaccard_fdr_a_vs_b" %in% colnames(t3_ok)) {
  png(file.path(plot_dir, "08_test3_fdr_comparison.png"),
      width = 1200, height = 800, res = 150)

  boxplot(
    list(
      "FDR+logFC\nA vs B" = t3_ok$jaccard_a_vs_b,
      "FDR only\nA vs B" = t3_ok$jaccard_fdr_a_vs_b,
      "FDR+logFC\nhalf vs full" = c(t3_ok$jaccard_a_vs_full,
                                     t3_ok$jaccard_b_vs_full),
      "FDR only\nhalf vs full" = c(t3_ok$jaccard_fdr_a_vs_full,
                                    t3_ok$jaccard_fdr_b_vs_full)
    ),
    col = c("lightyellow", "plum", "lightyellow", "plum"),
    main = "Split-half: FDR+logFC vs FDR-only",
    ylim = c(0, 1), ylab = "Jaccard"
  )

  dev.off()
  cat("08 saved\n")
}

# ============================================================
# 09: Summary — convergence curves
# ============================================================
if (has_t1 && has_t1b && has_t2) {
  png(file.path(plot_dir, "09_summary_convergence.png"),
      width = 2100, height = 700, res = 150)
  par(mfrow = c(1, 3), mar = c(5, 4, 3, 1))

  # Panel 1: Jaccard convergence with IQR
  t1_med <- aggregate(jaccard_vs_full ~ N_1st,
                       data = t1_ok, FUN = median)
  t1_q25 <- aggregate(jaccard_vs_full ~ N_1st, data = t1_ok,
                       FUN = function(x) quantile(x, 0.25))
  t1_q75 <- aggregate(jaccard_vs_full ~ N_1st, data = t1_ok,
                       FUN = function(x) quantile(x, 0.75))
  plot(t1_med$N_1st, t1_med$jaccard_vs_full,
       type = "b", pch = 19, ylim = c(0, 1),
       xlab = "N 1st-trim samples",
       ylab = "Jaccard vs full",
       main = "Jaccard convergence (IQR)")
  arrows(t1_med$N_1st, t1_q25$jaccard_vs_full,
         t1_med$N_1st, t1_q75$jaccard_vs_full,
         angle = 90, code = 3, length = 0.05, col = "grey40")

  # Panel 2: balanced vs unbalanced
  t2_med <- aggregate(jaccard_vs_full ~ N_total,
                       data = t2_ok, FUN = median)
  n_2nd <- t1_ok$N_total[1] - t1_ok$N_1st[1]
  plot(t1_med$N_1st + n_2nd, t1_med$jaccard_vs_full,
       type = "b", pch = 19, col = "blue", ylim = c(0, 1),
       xlab = "Total N", ylab = "Jaccard vs full",
       main = "Balanced (red) vs unbalanced (blue)")
  points(t2_med$N_total, t2_med$jaccard_vs_full,
         type = "b", pch = 17, col = "red")
  legend("bottomright",
         legend = c("Unbalanced (Test 1)",
                    "Balanced (Test 2)"),
         col = c("blue", "red"), pch = c(19, 17), bty = "n")

  # Panel 3: retention convergence
  t1_ret <- aggregate(overlap_vs_full ~ N_1st,
                       data = t1_ok, FUN = median)
  t1b_ret <- aggregate(overlap_vs_balanced ~ N_1st,
                        data = t1b_ok, FUN = median)
  plot(t1_ret$N_1st, t1_ret$overlap_vs_full,
       type = "b", pch = 19, col = "blue", ylim = c(0, 1),
       xlab = "N 1st-trim samples",
       ylab = "Median retention",
       main = "Retention convergence")
  points(t1b_ret$N_1st, t1b_ret$overlap_vs_balanced,
         type = "b", pch = 17, col = "red")
  legend("bottomright",
         legend = c(paste0("Full DEGs (", ref_n_deg, ")"),
                    paste0("Balanced DEGs (", bal_n_deg, ")")),
         col = c("blue", "red"), pch = c(19, 17), bty = "n")

  dev.off()
  cat("09 saved\n")
}

# ============================================================
# 10: Summary — full vs balanced comparison
# ============================================================
if (has_t1b) {
  png(file.path(plot_dir, "10_summary_full_vs_balanced.png"),
      width = 1800, height = 800, res = 150)
  par(mfrow = c(1, 2), mar = c(5, 4, 3, 1))

  # Panel 1: Jaccard vs full and vs balanced boxplots
  boxplot(jaccard_vs_full ~ N_1st, data = t1b_ok,
          col = "lightblue",
          xlab = "N first-trimester samples",
          ylab = "Jaccard index",
          main = "vs full integration", ylim = c(0, 1))

  boxplot(jaccard_vs_balanced ~ N_1st, data = t1b_ok,
          col = "lightsalmon",
          xlab = "N first-trimester samples",
          ylab = "Jaccard index",
          main = "vs balanced reference (2 datasets)",
          ylim = c(0, 1))

  dev.off()
  cat("10 saved\n")
}

cat("\nAll plots saved to:", plot_dir, "\n")
