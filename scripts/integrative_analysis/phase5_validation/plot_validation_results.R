#!/usr/bin/env Rscript
#
# Visualization for Phase 5 subsampling validation results.
# Reads TSVs produced by test1, test2, test3 scripts.
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

# ============================================================
# Test 1 plots
# ============================================================

t1_file <- file.path(input_dir, "test1_first_trim_subsample.tsv")
if (file.exists(t1_file)) {
  t1 <- read.delim(t1_file, stringsAsFactors = FALSE)
  t1_ok <- t1[t1$status == "ok", ]

  # DEG count
  png(file.path(plot_dir, "test1_deg_count.png"), width = 1200, height = 800, res = 150)
  boxplot(n_deg ~ N_1st, data = t1_ok, col = "lightblue",
          xlab = "N first-trimester samples", ylab = "Number of DEGs",
          main = "Test 1: DEG count vs first-trimester sample size")
  abline(h = ref_n_deg, lty = 2, col = "red", lwd = 2)
  legend("bottomright", legend = paste0("Full run (", ref_n_deg, " DEGs)"),
         lty = 2, col = "red", lwd = 2, bty = "n")
  dev.off()

  # Jaccard vs full
  png(file.path(plot_dir, "test1_jaccard_vs_full.png"), width = 1200, height = 800, res = 150)
  boxplot(jaccard_vs_full ~ N_1st, data = t1_ok, col = "lightgreen",
          xlab = "N first-trimester samples", ylab = "Jaccard index vs full run",
          main = "Test 1: DEG overlap with full result", ylim = c(0, 1))
  dev.off()

  # logFC correlation
  png(file.path(plot_dir, "test1_logfc_corr.png"), width = 1200, height = 800, res = 150)
  boxplot(logfc_pearson_vs_full ~ N_1st, data = t1_ok, col = "lightyellow",
          xlab = "N first-trimester samples", ylab = "Pearson r (logFC vs full)",
          main = "Test 1: Effect-size correlation with full result", ylim = c(0, 1))
  dev.off()

  # Within-size stability
  wj_file <- file.path(input_dir, "test1_within_size_jaccard.tsv")
  if (file.exists(wj_file)) {
    wj <- read.delim(wj_file, stringsAsFactors = FALSE)
    png(file.path(plot_dir, "test1_within_size_stability.png"),
        width = 1200, height = 800, res = 150)
    plot(wj$N_1st, wj$mean_pairwise_jaccard, type = "b", pch = 19,
         ylim = c(0, 1), xlab = "N first-trimester samples",
         ylab = "Mean pairwise Jaccard (within size)",
         main = "Test 1: Result stability across random draws")
    arrows(wj$N_1st,
           wj$mean_pairwise_jaccard - wj$sd_pairwise_jaccard,
           wj$N_1st,
           wj$mean_pairwise_jaccard + wj$sd_pairwise_jaccard,
           angle = 90, code = 3, length = 0.05, col = "grey40")
    dev.off()
  }
  # FDR-only: DEG count + Jaccard side by side
  if ("n_deg_fdr_only" %in% colnames(t1_ok)) {
    png(file.path(plot_dir, "test1_fdr_only_metrics.png"),
        width = 1800, height = 800, res = 150)
    par(mfrow = c(1, 3), mar = c(5, 4, 3, 1))

    boxplot(n_deg_fdr_only ~ N_1st, data = t1_ok, col = "lightblue",
            xlab = "N first-trimester samples",
            ylab = "Number of DEGs (FDR only)",
            main = "Test 1: FDR-only DEG count")
    abline(h = ref_n_deg_fdr, lty = 2, col = "red", lwd = 2)
    legend("bottomright",
           legend = paste0("Full (", ref_n_deg_fdr, " FDR-only)"),
           lty = 2, col = "red", lwd = 2, bty = "n")

    boxplot(jaccard_fdr_only ~ N_1st, data = t1_ok, col = "lightgreen",
            xlab = "N first-trimester samples",
            ylab = "Jaccard (FDR-only vs full FDR-only)",
            main = "Test 1: FDR-only overlap", ylim = c(0, 1))

    boxplot(jaccard_vs_full ~ N_1st, data = t1_ok, col = "lightyellow",
            xlab = "N first-trimester samples",
            ylab = "Jaccard vs full",
            main = "FDR+logFC for comparison", ylim = c(0, 1))

    dev.off()
  }
  cat("Test 1 plots saved\n")
}

# ============================================================
# Test 2 plots
# ============================================================

t2_file <- file.path(input_dir, "test2_balanced_subsample.tsv")
if (file.exists(t2_file)) {
  t2 <- read.delim(t2_file, stringsAsFactors = FALSE)
  t2_ok <- t2[t2$status == "ok", ]

  png(file.path(plot_dir, "test2_balanced_metrics.png"),
      width = 1800, height = 600, res = 150)
  par(mfrow = c(1, 3), mar = c(5, 4, 3, 1))

  boxplot(n_deg ~ n_per_trim, data = t2_ok, col = "lightsalmon",
          xlab = "N per trimester", ylab = "Number of DEGs",
          main = "Balanced: DEG count")
  abline(h = ref_n_deg, lty = 2, col = "red", lwd = 2)

  boxplot(jaccard_vs_full ~ n_per_trim, data = t2_ok, col = "lightsalmon",
          xlab = "N per trimester", ylab = "Jaccard vs full",
          main = "Balanced: overlap", ylim = c(0, 1))

  boxplot(logfc_pearson_vs_full ~ n_per_trim, data = t2_ok, col = "lightsalmon",
          xlab = "N per trimester", ylab = "Pearson r (logFC)",
          main = "Balanced: logFC correlation", ylim = c(0, 1))

  dev.off()

  if ("n_deg_fdr_only" %in% colnames(t2_ok)) {
    png(file.path(plot_dir, "test2_fdr_only_metrics.png"),
        width = 1200, height = 600, res = 150)
    par(mfrow = c(1, 2), mar = c(5, 4, 3, 1))

    boxplot(n_deg_fdr_only ~ n_per_trim, data = t2_ok,
            col = "lightsalmon",
            xlab = "N per trimester",
            ylab = "Number of DEGs (FDR only)",
            main = "Balanced: FDR-only DEG count")
    abline(h = ref_n_deg_fdr, lty = 2, col = "red", lwd = 2)

    boxplot(jaccard_fdr_only ~ n_per_trim, data = t2_ok,
            col = "lightsalmon",
            xlab = "N per trimester",
            ylab = "Jaccard (FDR-only vs full)",
            main = "Balanced: FDR-only overlap", ylim = c(0, 1))

    dev.off()
  }
  cat("Test 2 plots saved\n")
}

# ============================================================
# Test 3 plots
# ============================================================

t3_file <- file.path(input_dir, "test3_split_half.tsv")
if (file.exists(t3_file)) {
  t3 <- read.delim(t3_file, stringsAsFactors = FALSE)
  t3_ok <- t3[t3$status_a == "ok" & t3$status_b == "ok", ]

  png(file.path(plot_dir, "test3_split_half.png"),
      width = 1800, height = 600, res = 150)
  par(mfrow = c(1, 3), mar = c(5, 4, 3, 1))

  hist(t3_ok$jaccard_a_vs_b, breaks = 20, col = "plum",
       xlab = "Jaccard (half A vs half B)", main = "Split-half DEG overlap",
       xlim = c(0, 1))
  abline(v = median(t3_ok$jaccard_a_vs_b, na.rm = TRUE),
         lty = 2, col = "red", lwd = 2)

  hist(t3_ok$logfc_pearson_a_vs_b, breaks = 20, col = "plum",
       xlab = "Pearson r (half A vs half B)", main = "Split-half logFC correlation",
       xlim = c(0, 1))
  abline(v = median(t3_ok$logfc_pearson_a_vs_b, na.rm = TRUE),
         lty = 2, col = "red", lwd = 2)

  all_jaccard_vs_full <- c(t3_ok$jaccard_a_vs_full, t3_ok$jaccard_b_vs_full)
  hist(all_jaccard_vs_full, breaks = 20, col = "plum",
       xlab = "Jaccard (half vs full)", main = "Each half vs full run",
       xlim = c(0, 1))
  abline(v = median(all_jaccard_vs_full, na.rm = TRUE),
         lty = 2, col = "red", lwd = 2)

  dev.off()

  if ("jaccard_fdr_a_vs_b" %in% colnames(t3_ok)) {
    png(file.path(plot_dir, "test3_fdr_only_split_half.png"),
        width = 1800, height = 600, res = 150)
    par(mfrow = c(1, 3), mar = c(5, 4, 3, 1))

    hist(t3_ok$jaccard_fdr_a_vs_b, breaks = 20, col = "plum",
         xlab = "Jaccard FDR-only (A vs B)",
         main = "Split-half FDR-only overlap", xlim = c(0, 1))
    abline(v = median(t3_ok$jaccard_fdr_a_vs_b, na.rm = TRUE),
           lty = 2, col = "red", lwd = 2)

    fdr_vs_full <- c(t3_ok$jaccard_fdr_a_vs_full,
                     t3_ok$jaccard_fdr_b_vs_full)
    hist(fdr_vs_full, breaks = 20, col = "plum",
         xlab = "Jaccard FDR-only (half vs full)",
         main = "Each half vs full (FDR-only)", xlim = c(0, 1))
    abline(v = median(fdr_vs_full, na.rm = TRUE),
           lty = 2, col = "red", lwd = 2)

    boxplot(
      list(
        "FDR+logFC\nA vs B" = t3_ok$jaccard_a_vs_b,
        "FDR only\nA vs B" = t3_ok$jaccard_fdr_a_vs_b
      ),
      col = c("lightyellow", "plum"),
      main = "FDR+logFC vs FDR-only", ylim = c(0, 1),
      ylab = "Jaccard (A vs B)"
    )
    dev.off()
  }
  cat("Test 3 plots saved\n")
}

# ============================================================
# Combined summary
# ============================================================

has_all <- file.exists(t1_file) && file.exists(t2_file) && file.exists(t3_file)
if (has_all) {
  png(file.path(plot_dir, "combined_summary.png"),
      width = 2100, height = 700, res = 150)
  par(mfrow = c(1, 3), mar = c(5, 4, 3, 1))

  # Test 1: Jaccard convergence
  t1_med <- aggregate(jaccard_vs_full ~ N_1st, data = t1_ok, FUN = median)
  t1_q25 <- aggregate(jaccard_vs_full ~ N_1st, data = t1_ok,
                       FUN = function(x) quantile(x, 0.25, na.rm = TRUE))
  t1_q75 <- aggregate(jaccard_vs_full ~ N_1st, data = t1_ok,
                       FUN = function(x) quantile(x, 0.75, na.rm = TRUE))
  plot(t1_med$N_1st, t1_med$jaccard_vs_full, type = "b", pch = 19,
       ylim = c(0, 1), xlab = "N 1st-trim samples", ylab = "Jaccard vs full",
       main = "Test 1: Convergence")
  arrows(t1_med$N_1st, t1_q25$jaccard_vs_full, t1_med$N_1st, t1_q75$jaccard_vs_full,
         angle = 90, code = 3, length = 0.05, col = "grey40")

  # Test 2: Balanced vs unbalanced comparison
  t2_med <- aggregate(jaccard_vs_full ~ N_total, data = t2_ok, FUN = median)
  n_2nd <- t1_ok$N_total[1] - t1_ok$N_1st[1]
  plot(t1_med$N_1st + n_2nd, t1_med$jaccard_vs_full,
       type = "b", pch = 19, col = "blue", ylim = c(0, 1),
       xlab = "Total N", ylab = "Jaccard vs full",
       main = "Balanced (red) vs unbalanced (blue)")
  points(t2_med$N_total, t2_med$jaccard_vs_full, type = "b", pch = 17, col = "red")
  legend("bottomright", legend = c("Unbalanced (Test 1)", "Balanced (Test 2)"),
         col = c("blue", "red"), pch = c(19, 17), bty = "n")

  # Test 3: Split-half summary
  boxplot(
    list(
      "Jaccard\nA vs B" = t3_ok$jaccard_a_vs_b,
      "Jaccard\nhalf vs full" = all_jaccard_vs_full,
      "logFC r\nA vs B" = t3_ok$logfc_pearson_a_vs_b
    ),
    col = c("plum", "plum", "lightyellow"),
    main = "Test 3: Split-half metrics", ylim = c(0, 1)
  )

  dev.off()
  cat("Combined summary plot saved\n")

  # FDR-only vs FDR+logFC comparison across all tests
  has_fdr <- "jaccard_fdr_only" %in% colnames(t1_ok) &&
    "jaccard_fdr_only" %in% colnames(t2_ok) &&
    "jaccard_fdr_a_vs_b" %in% colnames(t3_ok)

  if (has_fdr) {
    png(file.path(plot_dir, "combined_fdr_comparison.png"),
        width = 2100, height = 700, res = 150)
    par(mfrow = c(1, 3), mar = c(5, 4, 3, 1))

    # Test 1: FDR-only vs FDR+logFC convergence
    t1_med_fdr <- aggregate(
      jaccard_fdr_only ~ N_1st, data = t1_ok, FUN = median
    )
    plot(t1_med$N_1st, t1_med$jaccard_vs_full,
         type = "b", pch = 19, col = "blue", ylim = c(0, 1),
         xlab = "N 1st-trim samples", ylab = "Jaccard vs full",
         main = "Test 1: FDR+logFC vs FDR-only")
    points(t1_med_fdr$N_1st, t1_med_fdr$jaccard_fdr_only,
           type = "b", pch = 17, col = "red")
    legend("bottomright",
           legend = c("FDR + logFC", "FDR only"),
           col = c("blue", "red"), pch = c(19, 17), bty = "n")

    # Test 2: Same comparison
    t2_med_fdr <- aggregate(
      jaccard_fdr_only ~ n_per_trim, data = t2_ok, FUN = median
    )
    plot(t2_med$N_total, t2_med$jaccard_vs_full,
         type = "b", pch = 19, col = "blue", ylim = c(0, 1),
         xlab = "Total N (balanced)", ylab = "Jaccard vs full",
         main = "Test 2: FDR+logFC vs FDR-only")
    points(t2_med_fdr$n_per_trim * 2,
           t2_med_fdr$jaccard_fdr_only,
           type = "b", pch = 17, col = "red")
    legend("bottomright",
           legend = c("FDR + logFC", "FDR only"),
           col = c("blue", "red"), pch = c(19, 17), bty = "n")

    # Test 3: Side-by-side boxplots
    boxplot(
      list(
        "FDR+logFC\nA vs B" = t3_ok$jaccard_a_vs_b,
        "FDR only\nA vs B" = t3_ok$jaccard_fdr_a_vs_b,
        "FDR+logFC\nvs full" = all_jaccard_vs_full,
        "FDR only\nvs full" = c(
          t3_ok$jaccard_fdr_a_vs_full,
          t3_ok$jaccard_fdr_b_vs_full
        )
      ),
      col = c("lightyellow", "plum",
              "lightyellow", "plum"),
      main = "Test 3: Split-half comparison",
      ylim = c(0, 1), ylab = "Jaccard"
    )
    dev.off()
    cat("FDR comparison plot saved\n")
  }
}

cat("\nAll plots saved to:", plot_dir, "\n")
