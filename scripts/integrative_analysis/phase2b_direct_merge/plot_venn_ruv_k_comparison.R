#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(VennDiagram)
  library(grid)
})

setwd("/home/shivers/a/r/igea-r/scripts/integrative-gene-expression-analysis")

fdr_cutoff <- 0.05

load_degs <- function(path) {
  if (!file.exists(path)) { cat(sprintf("MISSING: %s\n", path)); return(character(0)) }
  d <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   check.names = FALSE, quote = "")
  as.character(d$gene[!is.na(d$adj.P.Val) & d$adj.P.Val < fdr_cutoff])
}

comparisons <- list(
  list(
    id    = "1_2_all_datasets",
    label = "1st vs 2nd Trimester — All Datasets",
    combat = "output/phase2b_combat/phase2b_1_2_all_datasets/difexp_none_combat.tsv",
    ruv_k2 = "output/phase2b_ruv/phase2b_1_2_all_datasets_ruv/difexp_none_ruv.tsv",
    ruv_k3 = "output/phase2b_ruv/phase2b_1_2_all_datasets_ruv_k3/difexp_none_ruv.tsv",
    ruv_k4 = "output/phase2b_ruv/phase2b_1_2_all_datasets_ruv_k4/difexp_none_ruv.tsv"
  ),
  list(
    id    = "2_3_all_datasets",
    label = "2nd Trim vs Term — All Datasets",
    combat = "output/phase2b_combat/phase2b_2_3_all_datasets/difexp_none_combat.tsv",
    ruv_k2 = "output/phase2b_ruv/phase2b_2_3_all_datasets_ruv/difexp_none_ruv.tsv",
    ruv_k3 = "output/phase2b_ruv/phase2b_2_3_all_datasets_ruv_k3/difexp_none_ruv.tsv",
    ruv_k4 = "output/phase2b_ruv/phase2b_2_3_all_datasets_ruv_k4/difexp_none_ruv.tsv"
  )
)

out_dir <- "articles/imputation_article/misc"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (cmp in comparisons) {
  degs_combat <- load_degs(cmp$combat)
  degs_k2     <- load_degs(cmp$ruv_k2)
  degs_k3     <- load_degs(cmp$ruv_k3)
  degs_k4     <- load_degs(cmp$ruv_k4)

  cat(sprintf("\n%s:\n", cmp$label))
  cat(sprintf("  ComBat DEGs:   %d\n", length(degs_combat)))
  cat(sprintf("  RUV k=2 DEGs:  %d\n", length(degs_k2)))
  cat(sprintf("  RUV k=3 DEGs:  %d\n", length(degs_k3)))
  cat(sprintf("  RUV k=4 DEGs:  %d\n", length(degs_k4)))

  # Venn 1: ComBat vs RUV k=3
  out1 <- file.path(out_dir, sprintf("venn_ruv_k3_vs_combat_%s.pdf", cmp$id))
  v1 <- venn.diagram(
    x = list("ComBat" = degs_combat, "RUV (k=3)" = degs_k3),
    filename = NULL, fill = c("#4BACC6", "#F79646"), alpha = 0.5,
    cat.cex = 1.1, cat.fontface = "bold", cex = 1.4, fontface = "bold",
    main = paste(cmp$label, "— k=3"), main.cex = 1.2, main.fontface = "bold",
    sub = sprintf("DEGs at FDR < %.2f (no logFC cutoff)", fdr_cutoff),
    sub.cex = 0.9, margin = 0.1
  )
  pdf(out1, width = 6, height = 5); grid.draw(v1); dev.off()
  cat(sprintf("  Wrote: %s\n", out1))

  # Venn 2: ComBat vs RUV k=4
  out2 <- file.path(out_dir, sprintf("venn_ruv_k4_vs_combat_%s.pdf", cmp$id))
  v2 <- venn.diagram(
    x = list("ComBat" = degs_combat, "RUV (k=4)" = degs_k4),
    filename = NULL, fill = c("#4BACC6", "#C0504D"), alpha = 0.5,
    cat.cex = 1.1, cat.fontface = "bold", cex = 1.4, fontface = "bold",
    main = paste(cmp$label, "— k=4"), main.cex = 1.2, main.fontface = "bold",
    sub = sprintf("DEGs at FDR < %.2f (no logFC cutoff)", fdr_cutoff),
    sub.cex = 0.9, margin = 0.1
  )
  pdf(out2, width = 6, height = 5); grid.draw(v2); dev.off()
  cat(sprintf("  Wrote: %s\n", out2))

  # Venn 3: RUV k=2 vs k=3 vs k=4
  out3 <- file.path(out_dir, sprintf("venn_ruv_k2_k3_k4_%s.pdf", cmp$id))
  v3 <- venn.diagram(
    x = list("RUV k=2" = degs_k2, "RUV k=3" = degs_k3, "RUV k=4" = degs_k4),
    filename = NULL, fill = c("#9BBB59", "#F79646", "#C0504D"), alpha = 0.5,
    cat.cex = 1.0, cat.fontface = "bold", cex = 1.2, fontface = "bold",
    main = paste(cmp$label, "— RUV k comparison"), main.cex = 1.2, main.fontface = "bold",
    sub = sprintf("DEGs at FDR < %.2f (no logFC cutoff)", fdr_cutoff),
    sub.cex = 0.9, margin = 0.1
  )
  pdf(out3, width = 7, height = 6); grid.draw(v3); dev.off()
  cat(sprintf("  Wrote: %s\n", out3))
}
