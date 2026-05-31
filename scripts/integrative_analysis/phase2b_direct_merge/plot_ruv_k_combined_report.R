#!/usr/bin/env Rscript
#' Combined RUV k-comparison report: Venn diagrams + PCA for all runs.
#' Produces a single multi-page PDF with all diagnostics.

suppressPackageStartupMessages({
  library(VennDiagram)
  library(grid)
  library(ggplot2)
  library(patchwork)
})

setwd("/home/shivers/a/r/igea-r/scripts/integrative-gene-expression-analysis")

fdr_cutoff <- 0.05

pdata <- read.csv("data/phenodata/samples.csv", stringsAsFactors = FALSE)

load_degs <- function(path) {
  if (!file.exists(path)) { cat(sprintf("MISSING: %s\n", path)); return(character(0)) }
  d <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   check.names = FALSE, quote = "")
  as.character(d$gene[!is.na(d$adj.P.Val) & d$adj.P.Val < fdr_cutoff])
}

load_degs_with_logfc <- function(path, logfc_cutoff = 1) {
  if (!file.exists(path)) return(character(0))
  d <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   check.names = FALSE, quote = "")
  as.character(d$gene[!is.na(d$adj.P.Val) & d$adj.P.Val < fdr_cutoff & abs(d$logFC) > logfc_cutoff])
}

ruv_correct_matrix <- function(exprs, control_genes, k) {
  exprs <- as.matrix(exprs)
  control_idx <- which(rownames(exprs) %in% control_genes)
  if (length(control_idx) < k + 1) return(exprs)
  Y_c <- t(exprs[control_idx, , drop = FALSE])
  Y_c <- scale(Y_c, center = TRUE, scale = FALSE)
  svd_c <- svd(Y_c, nu = k, nv = 0)
  W <- svd_c$u[, seq_len(k), drop = FALSE]
  Yt <- t(exprs)
  alpha <- solve(t(W) %*% W) %*% t(W) %*% Yt
  t(Yt - W %*% alpha)
}

bruv_correct_matrix <- function(exprs, control_genes, k, sample_dataset, sample_group) {
  exprs <- as.matrix(exprs)
  control_idx <- which(rownames(exprs) %in% control_genes)
  if (length(control_idx) < k + 1) return(exprs)
  ds <- as.character(sample_dataset)
  grp <- as.character(sample_group)
  ds_groups <- split(grp, ds)
  balanced_ds <- names(which(sapply(ds_groups, function(g) length(unique(g)) > 1)))
  balanced_idx <- which(ds %in% balanced_ds)
  if (length(balanced_idx) < k + 2) return(exprs)
  Y_c_bal <- t(exprs[control_idx, balanced_idx, drop = FALSE])
  Y_c_bal <- scale(Y_c_bal, center = TRUE, scale = FALSE)
  svd_bal <- svd(Y_c_bal, nu = k, nv = k)
  bal_center <- colMeans(t(exprs[control_idx, balanced_idx, drop = FALSE]))
  Y_c_all <- sweep(t(exprs[control_idx, , drop = FALSE]), 2, bal_center, "-")
  V_k <- svd_bal$v[, seq_len(k), drop = FALSE]
  d_k <- svd_bal$d[seq_len(k)]
  W <- Y_c_all %*% V_k %*% diag(1 / d_k, nrow = k, ncol = k)
  Yt <- t(exprs)
  alpha <- solve(t(W) %*% W) %*% t(W) %*% Yt
  t(Yt - W %*% alpha)
}

make_pca_plot <- function(exprs, sample_meta, color_var, title, color_label = NULL) {
  exprs <- as.matrix(exprs)
  gene_vars <- apply(exprs, 1, var, na.rm = TRUE)
  keep <- gene_vars > 0 & !is.na(gene_vars)
  exprs <- exprs[keep, ]
  pca <- prcomp(t(exprs), center = TRUE, scale. = TRUE)
  var_pct <- summary(pca)$importance[2, 1:2] * 100
  df <- data.frame(
    PC1 = pca$x[, 1], PC2 = pca$x[, 2],
    color = as.factor(sample_meta[[color_var]][match(colnames(exprs), sample_meta$arraydatafile_exprscolumnnames)])
  )
  df <- df[!is.na(df$color), ]
  if (is.null(color_label)) color_label <- color_var
  ggplot(df, aes(x = PC1, y = PC2, color = color)) +
    geom_point(size = 2, alpha = 0.7) +
    stat_ellipse(type = "norm", level = 0.95, linetype = 2, show.legend = FALSE) +
    labs(title = title, x = sprintf("PC1 (%.1f%%)", var_pct[1]),
         y = sprintf("PC2 (%.1f%%)", var_pct[2]), color = color_label) +
    theme_bw(base_size = 9) +
    theme(plot.title = element_text(size = 9, face = "bold"),
          legend.position = "bottom",
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 7))
}

comparisons <- list(
  list(
    id = "1_2_balanced", label = "1st vs 2nd Trim — Balanced",
    combat_de = "output/phase2b_combat/phase2b_1_2_balanced/difexp_none_combat.tsv",
    ruv_de = list(
      k2 = "output/phase2b_ruv/phase2b_1_2_balanced_ruv/difexp_none_ruv.tsv",
      k3 = "output/phase2b_ruv/phase2b_1_2_balanced_ruv_k3/difexp_none_ruv.tsv",
      k4 = "output/phase2b_ruv/phase2b_1_2_balanced_ruv_k4/difexp_none_ruv.tsv",
      inv = "output/phase2b_ruv/phase2b_1_2_balanced_ruvinv/difexp_none_ruvinv.tsv",
      bruv = "output/phase2b_ruv/phase2b_1_2_balanced_bruv/difexp_none_bruv.tsv"
    ),
    uncorr = "output/phase2b_ruv/phase2b_1_2_balanced_ruv/exprs_none_ruv.tsv",
    combat_exprs = "output/phase2b_combat/phase2b_1_2_balanced/exprs_none_combat.tsv",
    ctrl_file = "articles/imputation_article/ruv_control_genes_1_2.txt",
    group_col = "Gestational.Age.Category",
    groups = c("First Trimester", "Second Trimester")
  ),
  list(
    id = "2_3_balanced", label = "2nd Trim vs Term — Balanced",
    combat_de = "output/phase2b_combat/phase2b_2_3_balanced/difexp_none_combat.tsv",
    ruv_de = list(
      k2 = "output/phase2b_ruv/phase2b_2_3_balanced_ruv/difexp_none_ruv.tsv",
      k3 = "output/phase2b_ruv/phase2b_2_3_balanced_ruv_k3/difexp_none_ruv.tsv",
      k4 = "output/phase2b_ruv/phase2b_2_3_balanced_ruv_k4/difexp_none_ruv.tsv",
      inv = "output/phase2b_ruv/phase2b_2_3_balanced_ruvinv/difexp_none_ruvinv.tsv",
      bruv = "output/phase2b_ruv/phase2b_2_3_balanced_bruv/difexp_none_bruv.tsv"
    ),
    uncorr = "output/phase2b_ruv/phase2b_2_3_balanced_ruv/exprs_none_ruv.tsv",
    combat_exprs = "output/phase2b_combat/phase2b_2_3_balanced/exprs_none_combat.tsv",
    ctrl_file = "articles/imputation_article/ruv_control_genes_2_3.txt",
    group_col = "Gestational.Age.Category",
    groups = c("Second Trimester", "Term")
  ),
  list(
    id = "1_2_all_datasets", label = "1st vs 2nd Trim — All Datasets",
    combat_de = "output/phase2b_combat/phase2b_1_2_all_datasets/difexp_none_combat.tsv",
    ruv_de = list(
      k2 = "output/phase2b_ruv/phase2b_1_2_all_datasets_ruv/difexp_none_ruv.tsv",
      k3 = "output/phase2b_ruv/phase2b_1_2_all_datasets_ruv_k3/difexp_none_ruv.tsv",
      k4 = "output/phase2b_ruv/phase2b_1_2_all_datasets_ruv_k4/difexp_none_ruv.tsv",
      inv = "output/phase2b_ruv/phase2b_1_2_all_datasets_ruvinv/difexp_none_ruvinv.tsv",
      bruv = "output/phase2b_ruv/phase2b_1_2_all_datasets_bruv/difexp_none_bruv.tsv"
    ),
    uncorr = "output/phase2b_ruv/phase2b_1_2_all_datasets_ruv/exprs_none_ruv.tsv",
    combat_exprs = "output/phase2b_combat/phase2b_1_2_all_datasets/exprs_none_combat.tsv",
    ctrl_file = "articles/imputation_article/ruv_control_genes_1_2.txt",
    group_col = "Gestational.Age.Category",
    groups = c("First Trimester", "Second Trimester")
  ),
  list(
    id = "2_3_all_datasets", label = "2nd Trim vs Term — All Datasets",
    combat_de = "output/phase2b_combat/phase2b_2_3_all_datasets/difexp_none_combat.tsv",
    ruv_de = list(
      k2 = "output/phase2b_ruv/phase2b_2_3_all_datasets_ruv/difexp_none_ruv.tsv",
      k3 = "output/phase2b_ruv/phase2b_2_3_all_datasets_ruv_k3/difexp_none_ruv.tsv",
      k4 = "output/phase2b_ruv/phase2b_2_3_all_datasets_ruv_k4/difexp_none_ruv.tsv",
      inv = "output/phase2b_ruv/phase2b_2_3_all_datasets_ruvinv/difexp_none_ruvinv.tsv",
      bruv = "output/phase2b_ruv/phase2b_2_3_all_datasets_bruv/difexp_none_bruv.tsv"
    ),
    uncorr = "output/phase2b_ruv/phase2b_2_3_all_datasets_ruv/exprs_none_ruv.tsv",
    combat_exprs = "output/phase2b_combat/phase2b_2_3_all_datasets/exprs_none_combat.tsv",
    ctrl_file = "articles/imputation_article/ruv_control_genes_2_3.txt",
    group_col = "Gestational.Age.Category",
    groups = c("Second Trimester", "Term")
  )
)

out_dir <- "articles/imputation_article/misc"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_pdf <- file.path(out_dir, "ruv_k_comparison_report.pdf")

pdf(out_pdf, width = 16, height = 10)

for (cmp in comparisons) {
  cat(sprintf("\n=== %s ===\n", cmp$label))

  # --- Venn diagrams page ---
  degs_combat <- load_degs(cmp$combat_de)
  degs_k2 <- load_degs(cmp$ruv_de$k2)
  degs_k3 <- load_degs(cmp$ruv_de$k3)
  degs_k4 <- load_degs(cmp$ruv_de$k4)
  degs_inv <- load_degs(cmp$ruv_de$inv)

  degs_combat_lfc <- load_degs_with_logfc(cmp$combat_de)
  degs_k2_lfc <- load_degs_with_logfc(cmp$ruv_de$k2)
  degs_k3_lfc <- load_degs_with_logfc(cmp$ruv_de$k3)
  degs_k4_lfc <- load_degs_with_logfc(cmp$ruv_de$k4)
  degs_inv_lfc <- load_degs_with_logfc(cmp$ruv_de$inv)
  degs_bruv <- load_degs(cmp$ruv_de$bruv)
  degs_bruv_lfc <- load_degs_with_logfc(cmp$ruv_de$bruv)

  cat(sprintf("  ComBat: %d (FDR) / %d (FDR+logFC)\n", length(degs_combat), length(degs_combat_lfc)))
  cat(sprintf("  RUV k=2: %d / %d\n", length(degs_k2), length(degs_k2_lfc)))
  cat(sprintf("  RUV k=3: %d / %d\n", length(degs_k3), length(degs_k3_lfc)))
  cat(sprintf("  RUV k=4: %d / %d\n", length(degs_k4), length(degs_k4_lfc)))
  cat(sprintf("  RUVinv:  %d / %d\n", length(degs_inv), length(degs_inv_lfc)))
  cat(sprintf("  BRUV:    %d / %d\n", length(degs_bruv), length(degs_bruv_lfc)))

  # Page 1: Venn diagrams — ComBat vs each RUV variant
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(3, 2,
    heights = unit(c(0.08, 0.46, 0.46), "npc"))))

  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1:2))
  grid.text(paste(cmp$label, "— DEG Venn Diagrams (FDR < 0.05)"),
            gp = gpar(fontsize = 16, fontface = "bold"))
  popViewport()

  v1 <- venn.diagram(
    x = list("ComBat" = degs_combat, "RUV k=2" = degs_k2),
    filename = NULL, fill = c("#4BACC6", "#9BBB59"), alpha = 0.5,
    cat.cex = 1.0, cat.fontface = "bold", cex = 1.2, fontface = "bold",
    main = "ComBat vs RUV k=2", main.cex = 1.1, main.fontface = "bold",
    margin = 0.15)
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
  grid.draw(v1)
  popViewport()

  v2 <- venn.diagram(
    x = list("ComBat" = degs_combat, "RUV k=3" = degs_k3),
    filename = NULL, fill = c("#4BACC6", "#F79646"), alpha = 0.5,
    cat.cex = 1.0, cat.fontface = "bold", cex = 1.2, fontface = "bold",
    main = "ComBat vs RUV k=3", main.cex = 1.1, main.fontface = "bold",
    margin = 0.15)
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
  grid.draw(v2)
  popViewport()

  v3 <- venn.diagram(
    x = list("ComBat" = degs_combat, "RUV k=4" = degs_k4),
    filename = NULL, fill = c("#4BACC6", "#C0504D"), alpha = 0.5,
    cat.cex = 1.0, cat.fontface = "bold", cex = 1.2, fontface = "bold",
    main = "ComBat vs RUV k=4", main.cex = 1.1, main.fontface = "bold",
    margin = 0.15)
  pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1))
  grid.draw(v3)
  popViewport()

  v4 <- venn.diagram(
    x = list("ComBat" = degs_combat, "RUVinv" = degs_inv),
    filename = NULL, fill = c("#4BACC6", "#7030A0"), alpha = 0.5,
    cat.cex = 1.0, cat.fontface = "bold", cex = 1.2, fontface = "bold",
    main = "ComBat vs RUVinv", main.cex = 1.1, main.fontface = "bold",
    margin = 0.15)
  pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 2))
  grid.draw(v4)
  popViewport()
  popViewport()

  # Page 1b: BRUV Venns
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(3, 2,
    heights = unit(c(0.08, 0.46, 0.46), "npc"))))

  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1:2))
  grid.text(paste(cmp$label, "— BRUV Venn Diagrams (FDR < 0.05)"),
            gp = gpar(fontsize = 16, fontface = "bold"))
  popViewport()

  vb1 <- venn.diagram(
    x = list("ComBat" = degs_combat, "BRUV" = degs_bruv),
    filename = NULL, fill = c("#4BACC6", "#E07000"), alpha = 0.5,
    cat.cex = 1.0, cat.fontface = "bold", cex = 1.2, fontface = "bold",
    main = "ComBat vs BRUV", main.cex = 1.1, main.fontface = "bold",
    margin = 0.15)
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
  grid.draw(vb1)
  popViewport()

  vb2 <- venn.diagram(
    x = list("RUV k=2" = degs_k2, "BRUV" = degs_bruv),
    filename = NULL, fill = c("#9BBB59", "#E07000"), alpha = 0.5,
    cat.cex = 1.0, cat.fontface = "bold", cex = 1.2, fontface = "bold",
    main = "RUV k=2 vs BRUV", main.cex = 1.1, main.fontface = "bold",
    margin = 0.15)
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
  grid.draw(vb2)
  popViewport()

  vb3 <- venn.diagram(
    x = list("ComBat" = degs_combat, "RUV k=2" = degs_k2, "BRUV" = degs_bruv),
    filename = NULL, fill = c("#4BACC6", "#9BBB59", "#E07000"), alpha = 0.5,
    cat.cex = 1.0, cat.fontface = "bold", cex = 1.0, fontface = "bold",
    main = "ComBat vs RUV k=2 vs BRUV", main.cex = 1.1, main.fontface = "bold",
    margin = 0.15)
  pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1))
  grid.draw(vb3)
  popViewport()

  vb4 <- venn.diagram(
    x = list("ComBat" = degs_combat, "BRUV" = degs_bruv, "RUVinv" = degs_inv),
    filename = NULL, fill = c("#4BACC6", "#E07000", "#7030A0"), alpha = 0.5,
    cat.cex = 1.0, cat.fontface = "bold", cex = 1.0, fontface = "bold",
    main = "ComBat vs BRUV vs RUVinv", main.cex = 1.1, main.fontface = "bold",
    margin = 0.15)
  pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 2))
  grid.draw(vb4)
  popViewport()
  popViewport()

  # Page 2: Venn — RUV variants comparison
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(3, 2,
    heights = unit(c(0.08, 0.46, 0.46), "npc"))))

  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1:2))
  grid.text(paste(cmp$label, "— RUV Variant Comparison (FDR < 0.05)"),
            gp = gpar(fontsize = 16, fontface = "bold"))
  popViewport()

  v5 <- venn.diagram(
    x = list("k=2" = degs_k2, "k=3" = degs_k3, "k=4" = degs_k4),
    filename = NULL, fill = c("#9BBB59", "#F79646", "#C0504D"), alpha = 0.5,
    cat.cex = 1.0, cat.fontface = "bold", cex = 1.0, fontface = "bold",
    main = "RUV k=2 vs k=3 vs k=4", main.cex = 1.1, main.fontface = "bold",
    margin = 0.15)
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
  grid.draw(v5)
  popViewport()

  v6 <- venn.diagram(
    x = list("RUV k=2" = degs_k2, "RUVinv" = degs_inv),
    filename = NULL, fill = c("#9BBB59", "#7030A0"), alpha = 0.5,
    cat.cex = 1.0, cat.fontface = "bold", cex = 1.2, fontface = "bold",
    main = "RUV k=2 vs RUVinv", main.cex = 1.1, main.fontface = "bold",
    margin = 0.15)
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
  grid.draw(v6)
  popViewport()

  v7 <- venn.diagram(
    x = list("ComBat" = degs_combat, "RUV k=2" = degs_k2, "RUVinv" = degs_inv),
    filename = NULL, fill = c("#4BACC6", "#9BBB59", "#7030A0"), alpha = 0.5,
    cat.cex = 1.0, cat.fontface = "bold", cex = 1.0, fontface = "bold",
    main = "ComBat vs RUV k=2 vs RUVinv", main.cex = 1.1, main.fontface = "bold",
    margin = 0.15)
  pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1))
  grid.draw(v7)
  popViewport()
  popViewport()

  # --- PCA page ---
  if (!file.exists(cmp$uncorr) || !file.exists(cmp$combat_exprs)) {
    cat("  Skipping PCA: missing expression files\n")
    next
  }

  uncorr_exprs <- read.table(cmp$uncorr, header = TRUE, sep = "\t",
                              check.names = FALSE, row.names = 1)
  combat_exprs <- read.table(cmp$combat_exprs, header = TRUE, sep = "\t",
                              check.names = FALSE, row.names = 1)
  ctrl_genes <- readLines(cmp$ctrl_file)

  ruv_k2 <- ruv_correct_matrix(uncorr_exprs, ctrl_genes, 2)
  ruv_k3 <- ruv_correct_matrix(uncorr_exprs, ctrl_genes, 3)
  ruv_k4 <- ruv_correct_matrix(uncorr_exprs, ctrl_genes, 4)
  # RUVinv: use high k (rank of control gene space) for visualization
  n_ctrl <- sum(rownames(uncorr_exprs) %in% ctrl_genes)
  ruvinv_vis_k <- min(n_ctrl - 1, ncol(uncorr_exprs) - 1, 20)
  ruv_inv <- ruv_correct_matrix(uncorr_exprs, ctrl_genes, ruvinv_vis_k)

  samples <- colnames(uncorr_exprs)
  meta <- pdata[pdata$arraydatafile_exprscolumnnames %in% samples, ]
  meta <- meta[meta[[cmp$group_col]] %in% cmp$groups, ]
  common_genes <- intersect(rownames(uncorr_exprs), rownames(combat_exprs))

  p1 <- make_pca_plot(uncorr_exprs, meta, "secondaryaccession", "Uncorrected", "Dataset")
  p2 <- make_pca_plot(combat_exprs[common_genes, ], meta, "secondaryaccession", "ComBat", "Dataset")
  p3 <- make_pca_plot(ruv_k2, meta, "secondaryaccession", "RUV k=2", "Dataset")
  p4 <- make_pca_plot(ruv_k3, meta, "secondaryaccession", "RUV k=3", "Dataset")
  p5 <- make_pca_plot(ruv_k4, meta, "secondaryaccession", "RUV k=4", "Dataset")

  p6  <- make_pca_plot(uncorr_exprs, meta, cmp$group_col, "Uncorrected", "Trimester")
  p7  <- make_pca_plot(combat_exprs[common_genes, ], meta, cmp$group_col, "ComBat", "Trimester")
  p8  <- make_pca_plot(ruv_k2, meta, cmp$group_col, "RUV k=2", "Trimester")
  p9  <- make_pca_plot(ruv_k3, meta, cmp$group_col, "RUV k=3", "Trimester")
  p10 <- make_pca_plot(ruv_k4, meta, cmp$group_col, "RUV k=4", "Trimester")

  combined <- (p1 | p2 | p3 | p4 | p5) / (p6 | p7 | p8 | p9 | p10) +
    plot_annotation(
      title = paste(cmp$label, "— PCA (k=2,3,4)"),
      subtitle = "Top: colored by dataset (batch). Bottom: colored by trimester (biology).",
      theme = theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10)
      )
    )
  print(combined)

  # Extra PCA page: RUVinv vs ComBat
  p_inv1 <- make_pca_plot(ruv_inv, meta, "secondaryaccession",
                           sprintf("RUVinv (k=%d proj.)", ruvinv_vis_k), "Dataset")
  p_inv2 <- make_pca_plot(ruv_inv, meta, cmp$group_col,
                           sprintf("RUVinv (k=%d proj.)", ruvinv_vis_k), "Trimester")

  combined_inv <- (p1 | p2 | p3 | p_inv1) / (p6 | p7 | p8 | p_inv2) +
    plot_annotation(
      title = paste(cmp$label, "— PCA: Uncorrected / ComBat / RUV k=2 / RUVinv"),
      subtitle = "Top: dataset. Bottom: trimester. RUVinv visualized by projecting out top-20 SVD components.",
      theme = theme(
        plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 9)
      )
    )
  print(combined_inv)

  # BRUV PCA page
  sample_ds <- meta$secondaryaccession[match(colnames(uncorr_exprs), meta$arraydatafile_exprscolumnnames)]
  sample_grp <- meta[[cmp$group_col]][match(colnames(uncorr_exprs), meta$arraydatafile_exprscolumnnames)]
  bruv_corrected <- bruv_correct_matrix(uncorr_exprs, ctrl_genes, 2, sample_ds, sample_grp)

  p_bruv1 <- make_pca_plot(bruv_corrected, meta, "secondaryaccession", "BRUV k=2", "Dataset")
  p_bruv2 <- make_pca_plot(bruv_corrected, meta, cmp$group_col, "BRUV k=2", "Trimester")

  combined_bruv <- (p1 | p2 | p3 | p_bruv1) / (p6 | p7 | p8 | p_bruv2) +
    plot_annotation(
      title = paste(cmp$label, "— PCA: Uncorrected / ComBat / RUV k=2 / BRUV"),
      subtitle = "Top: dataset. Bottom: trimester. BRUV estimates W from balanced datasets only, projects onto all.",
      theme = theme(
        plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 9)
      )
    )
  print(combined_bruv)
  cat(sprintf("  PCA pages added\n"))
}

dev.off()
cat(sprintf("\nWrote: %s\n", out_pdf))
