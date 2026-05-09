#!/usr/bin/env Rscript
#' Validate RUV batch correction via PCA.
#'
#' For each contrast (1_2, 2_3) x dataset selection (balanced, all_datasets),
#' produces a 2x3 panel of PCA plots:
#'   row 1: colored by dataset (batch)
#'   row 2: colored by trimester (biology)
#'   columns: Uncorrected | ComBat | RUV-corrected
#'
#' The RUV-corrected matrix is computed as Y - W %*% alpha (projection of
#' unwanted factors removed from expression) for visualization only.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

setwd("/home/shivers/a/r/igea-r/scripts/integrative-gene-expression-analysis")

runs <- list(
  list(
    id       = "1_2_balanced",
    label    = "1st vs 2nd Trimester - Balanced",
    uncorr   = "output/phase2b_ruv/phase2b_1_2_balanced_ruv/exprs_none_ruv.tsv",
    combat   = "output/phase2b_combat/phase2b_1_2_balanced/exprs_none_combat.tsv",
    ruv_exprs = "output/phase2b_ruv/phase2b_1_2_balanced_ruv/exprs_none_ruv.tsv",
    ruv_ctrl  = "articles/imputation_article/ruv_control_genes_1_2.txt",
    ruv_k     = 2,
    group_col = "Gestational.Age.Category",
    groups    = c("First Trimester", "Second Trimester")
  ),
  list(
    id       = "2_3_balanced",
    label    = "2nd Trim vs Term - Balanced",
    uncorr   = "output/phase2b_ruv/phase2b_2_3_balanced_ruv/exprs_none_ruv.tsv",
    combat   = "output/phase2b_combat/phase2b_2_3_balanced/exprs_none_combat.tsv",
    ruv_exprs = "output/phase2b_ruv/phase2b_2_3_balanced_ruv/exprs_none_ruv.tsv",
    ruv_ctrl  = "articles/imputation_article/ruv_control_genes_2_3.txt",
    ruv_k     = 2,
    group_col = "Gestational.Age.Category",
    groups    = c("Second Trimester", "Term")
  ),
  list(
    id       = "1_2_all_datasets",
    label    = "1st vs 2nd Trimester - All Datasets",
    uncorr   = "output/phase2b_ruv/phase2b_1_2_all_datasets_ruv/exprs_none_ruv.tsv",
    combat   = "output/phase2b_combat/phase2b_1_2_all_datasets/exprs_none_combat.tsv",
    ruv_exprs = "output/phase2b_ruv/phase2b_1_2_all_datasets_ruv/exprs_none_ruv.tsv",
    ruv_ctrl  = "articles/imputation_article/ruv_control_genes_1_2.txt",
    ruv_k     = 2,
    group_col = "Gestational.Age.Category",
    groups    = c("First Trimester", "Second Trimester")
  ),
  list(
    id       = "2_3_all_datasets",
    label    = "2nd Trim vs Term - All Datasets",
    uncorr   = "output/phase2b_ruv/phase2b_2_3_all_datasets_ruv/exprs_none_ruv.tsv",
    combat   = "output/phase2b_combat/phase2b_2_3_all_datasets/exprs_none_combat.tsv",
    ruv_exprs = "output/phase2b_ruv/phase2b_2_3_all_datasets_ruv/exprs_none_ruv.tsv",
    ruv_ctrl  = "articles/imputation_article/ruv_control_genes_2_3.txt",
    ruv_k     = 2,
    group_col = "Gestational.Age.Category",
    groups    = c("Second Trimester", "Term")
  )
)

pdata <- read.csv("data/phenodata/samples.csv", stringsAsFactors = FALSE)

ruv_correct_matrix <- function(exprs, control_genes, k) {
  exprs <- as.matrix(exprs)
  control_idx <- which(rownames(exprs) %in% control_genes)
  if (length(control_idx) < k + 1) {
    cat(sprintf("  Warning: only %d control genes, need %d\n", length(control_idx), k + 1))
    return(exprs)
  }
  Y_c <- t(exprs[control_idx, , drop = FALSE])
  Y_c <- scale(Y_c, center = TRUE, scale = FALSE)
  svd_c <- svd(Y_c, nu = k, nv = 0)
  W <- svd_c$u[, seq_len(k), drop = FALSE]

  # Project out W from the full matrix: Y_corr = Y - W %*% (W+ %*% Y)
  # where W+ = (W'W)^{-1} W' is the pseudoinverse
  Yt <- t(exprs)  # samples x genes
  alpha <- solve(t(W) %*% W) %*% t(W) %*% Yt
  corrected <- t(Yt - W %*% alpha)
  corrected
}

make_pca_plot <- function(exprs, sample_meta, color_var, title, color_label = NULL) {
  exprs <- as.matrix(exprs)
  gene_vars <- apply(exprs, 1, var, na.rm = TRUE)
  keep <- gene_vars > 0 & !is.na(gene_vars)
  exprs <- exprs[keep, ]

  pca <- prcomp(t(exprs), center = TRUE, scale. = TRUE)
  var_pct <- summary(pca)$importance[2, 1:2] * 100

  df <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    color = as.factor(sample_meta[[color_var]][match(colnames(exprs), sample_meta$arraydatafile_exprscolumnnames)])
  )
  df <- df[!is.na(df$color), ]

  if (is.null(color_label)) color_label <- color_var

  ggplot(df, aes(x = PC1, y = PC2, color = color)) +
    geom_point(size = 2, alpha = 0.7) +
    stat_ellipse(type = "norm", level = 0.95, linetype = 2, show.legend = FALSE) +
    labs(
      title = title,
      x = sprintf("PC1 (%.1f%%)", var_pct[1]),
      y = sprintf("PC2 (%.1f%%)", var_pct[2]),
      color = color_label
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      legend.position = "bottom",
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 8)
    )
}

out_dir <- "articles/imputation_article/misc"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (run in runs) {
  cat(sprintf("\n=== %s ===\n", run$label))

  if (!file.exists(run$uncorr) || !file.exists(run$combat)) {
    cat("  Skipping: missing files\n")
    next
  }

  uncorr_exprs <- read.table(run$uncorr, header = TRUE, sep = "\t",
                              check.names = FALSE, row.names = 1)
  combat_exprs <- read.table(run$combat, header = TRUE, sep = "\t",
                              check.names = FALSE, row.names = 1)

  ctrl_genes <- readLines(run$ruv_ctrl)
  ruv_exprs <- ruv_correct_matrix(uncorr_exprs, ctrl_genes, run$ruv_k)

  # Subset pdata to samples in the matrices
  samples <- colnames(uncorr_exprs)
  meta <- pdata[pdata$arraydatafile_exprscolumnnames %in% samples, ]
  meta <- meta[meta[[run$group_col]] %in% run$groups, ]
  cat(sprintf("  Samples: %d, Genes: %d\n", length(samples), nrow(uncorr_exprs)))

  # Common genes across uncorrected and combat (combat may drop zero-var genes)
  common_genes <- intersect(rownames(uncorr_exprs), rownames(combat_exprs))
  cat(sprintf("  Common genes (uncorr & combat): %d\n", length(common_genes)))

  # Row 1: colored by dataset
  p1 <- make_pca_plot(uncorr_exprs, meta, "secondaryaccession", "Uncorrected", "Dataset")
  p2 <- make_pca_plot(combat_exprs[common_genes, ], meta, "secondaryaccession", "ComBat", "Dataset")
  p3 <- make_pca_plot(ruv_exprs, meta, "secondaryaccession", "RUV (k=2)", "Dataset")

  # Row 2: colored by trimester
  p4 <- make_pca_plot(uncorr_exprs, meta, run$group_col, "Uncorrected", "Trimester")
  p5 <- make_pca_plot(combat_exprs[common_genes, ], meta, run$group_col, "ComBat", "Trimester")
  p6 <- make_pca_plot(ruv_exprs, meta, run$group_col, "RUV (k=2)", "Trimester")

  combined <- (p1 | p2 | p3) / (p4 | p5 | p6) +
    plot_annotation(
      title = run$label,
      subtitle = "Top row: colored by dataset (batch). Bottom row: colored by trimester (biology).",
      theme = theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10)
      )
    )

  out_file <- file.path(out_dir, sprintf("pca_ruv_validation_%s.pdf", run$id))
  ggsave(out_file, combined, width = 16, height = 10)
  cat(sprintf("  Wrote: %s\n", out_file))

  out_png <- file.path(out_dir, sprintf("pca_ruv_validation_%s.png", run$id))
  png(out_png, width = 16, height = 10, units = "in", res = 150, type = "cairo")
  print(combined)
  dev.off()
  cat(sprintf("  Wrote: %s\n", out_png))
}

cat("\nDone.\n")
