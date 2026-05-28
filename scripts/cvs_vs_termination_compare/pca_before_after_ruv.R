#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

setwd("/home/shivers/a/r/igea-r/scripts/integrative-gene-expression-analysis")

pheno <- read.csv("data/phenodata/samples_cvs_ga_matched.csv",
                   stringsAsFactors = FALSE, quote = "\"")
cvs_expr <- read.delim("data/mapped/cvs/GSE12767_entrez_protein_coding.tsv",
                        row.names = 1, check.names = FALSE)
cvs_pheno <- pheno[pheno$secondaryaccession == "GSE12767" &
                    pheno$Diagnosis == "Healthy", ]
cvs_expr <- cvs_expr[, colnames(cvs_expr) %in% cvs_pheno$arraydatafile_exprscolumnnames]

config_dir <- "scripts/integrative_analysis/phase2b_direct_merge/config_cvs_vs_termination_compare"

comparisons <- list(
  list(name = "GSE100051", file = "data/mapped/GSE100051.tsv",
       ga = c("10", "11", "12"), platform = "Illumina HT-12"),
  list(name = "GSE93520",  file = "data/mapped/GSE93520.tsv",
       ga = c("10"), platform = "Agilent 4x44K"),
  list(name = "GSE28551",  file = "data/mapped/GSE28551.tsv",
       ga = NULL, platform = "ABI Human Genome v2")
)

for (comp in comparisons) {
  cat(sprintf("\n=== PCA: CVS vs %s ===\n", comp$name))

  abort_expr <- read.delim(comp$file, row.names = 1, check.names = FALSE)

  if (!is.null(comp$ga)) {
    ap <- pheno[pheno$secondaryaccession == comp$name &
                pheno$Diagnosis == "Healthy" &
                pheno$Gestational.Age %in% comp$ga, ]
  } else {
    ap <- pheno[pheno$secondaryaccession == comp$name &
                pheno$Diagnosis == "Healthy" &
                pheno$Gestational.Age.Category == "First Trimester", ]
  }
  abort_expr <- abort_expr[, colnames(abort_expr) %in% ap$arraydatafile_exprscolumnnames]

  shared <- intersect(rownames(cvs_expr), rownames(abort_expr))
  merged <- cbind(cvs_expr[shared, ], abort_expr[shared, ])

  gene_types <- suppressMessages(
    AnnotationDbi::select(org.Hs.eg.db, keys = rownames(merged),
                          columns = "GENETYPE", keytype = "ENTREZID")
  )
  pc <- gene_types$ENTREZID[gene_types$GENETYPE == "protein-coding" & !is.na(gene_types$GENETYPE)]
  merged <- merged[rownames(merged) %in% pc, ]
  vars <- apply(merged, 1, var)
  merged <- merged[vars > 0.01, ]

  n_cvs <- ncol(cvs_expr)
  n_abort <- ncol(abort_expr)

  sample_info <- data.frame(
    sample = colnames(merged),
    dataset = c(rep("GSE12767", n_cvs), rep(comp$name, n_abort)),
    group = c(rep("CVS", n_cvs), rep("Termination", n_abort)),
    stringsAsFactors = FALSE
  )

  control_genes <- readLines(file.path(config_dir,
    sprintf("control_genes_cvs_%s.txt", comp$name)))
  control_idx <- which(rownames(merged) %in% control_genes)

  # ── Before RUV PCA ──
  pca_before <- prcomp(t(merged), center = TRUE, scale. = FALSE)
  var_before <- round(100 * pca_before$sdev^2 / sum(pca_before$sdev^2), 1)

  df_before <- data.frame(
    PC1 = pca_before$x[, 1], PC2 = pca_before$x[, 2],
    Dataset = sample_info$dataset, Group = sample_info$group
  )

  # ── RUV correction for visualization ──
  merged_mat <- as.matrix(merged)
  Y_c <- t(merged_mat[control_idx, , drop = FALSE])
  Y_c <- scale(Y_c, center = TRUE, scale = FALSE)
  svd_c <- svd(Y_c, nu = 1, nv = 0)
  W <- svd_c$u[, 1, drop = FALSE]
  alpha <- solve(t(W) %*% W) %*% t(W) %*% t(merged_mat)
  corrected <- merged_mat - t(W %*% alpha)

  pca_after <- prcomp(t(corrected), center = TRUE, scale. = FALSE)
  var_after <- round(100 * pca_after$sdev^2 / sum(pca_after$sdev^2), 1)

  df_after <- data.frame(
    PC1 = pca_after$x[, 1], PC2 = pca_after$x[, 2],
    Dataset = sample_info$dataset, Group = sample_info$group
  )

  # ── Plot ──
  out_dir <- sprintf("output/cvs_vs_termination_compare/cvs_vs_%s", comp$name)
  png_file <- file.path(out_dir, "pca_before_after_ruv.png")

  png(png_file, width = 14, height = 6, units = "in", res = 200)
  par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))

  # Before
  cols_before <- ifelse(df_before$Dataset == "GSE12767", "#1565C0", "#E53935")
  pchs_before <- ifelse(df_before$Group == "CVS", 16, 17)
  plot(df_before$PC1, df_before$PC2,
       col = cols_before, pch = pchs_before, cex = 1.5,
       xlab = sprintf("PC1 (%.1f%%)", var_before[1]),
       ylab = sprintf("PC2 (%.1f%%)", var_before[2]),
       main = sprintf("Before RUV — CVS vs %s", comp$name))
  legend("topright",
         legend = c("GSE12767 (CVS)", sprintf("%s (Term.)", comp$name)),
         col = c("#1565C0", "#E53935"), pch = c(16, 17), cex = 0.8)

  # After
  cols_after <- ifelse(df_after$Dataset == "GSE12767", "#1565C0", "#E53935")
  pchs_after <- ifelse(df_after$Group == "CVS", 16, 17)
  plot(df_after$PC1, df_after$PC2,
       col = cols_after, pch = pchs_after, cex = 1.5,
       xlab = sprintf("PC1 (%.1f%%)", var_after[1]),
       ylab = sprintf("PC2 (%.1f%%)", var_after[2]),
       main = sprintf("After RUV k=1 — CVS vs %s", comp$name))
  legend("topright",
         legend = c("GSE12767 (CVS)", sprintf("%s (Term.)", comp$name)),
         col = c("#1565C0", "#E53935"), pch = c(16, 17), cex = 0.8)

  dev.off()
  cat(sprintf("  Saved: %s\n", png_file))
  cat(sprintf("  Before: PC1=%.1f%%, PC2=%.1f%%\n", var_before[1], var_before[2]))
  cat(sprintf("  After:  PC1=%.1f%%, PC2=%.1f%%\n", var_after[1], var_after[2]))
}

cat("\nPCA plots complete.\n")
