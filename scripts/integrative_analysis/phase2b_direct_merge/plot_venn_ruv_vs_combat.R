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
    id    = "1_2_balanced",
    label = "1st vs 2nd Trimester — Balanced",
    combat = "output/phase2b_combat/phase2b_1_2_balanced/difexp_none_combat.tsv",
    ruv    = "output/phase2b_ruv/phase2b_1_2_balanced_ruv/difexp_none_ruv.tsv"
  ),
  list(
    id    = "2_3_balanced",
    label = "2nd Trim vs Term — Balanced",
    combat = "output/phase2b_combat/phase2b_2_3_balanced/difexp_none_combat.tsv",
    ruv    = "output/phase2b_ruv/phase2b_2_3_balanced_ruv/difexp_none_ruv.tsv"
  ),
  list(
    id    = "1_2_all_datasets",
    label = "1st vs 2nd Trimester — All Datasets",
    combat = "output/phase2b_combat/phase2b_1_2_all_datasets/difexp_none_combat.tsv",
    ruv    = "output/phase2b_ruv/phase2b_1_2_all_datasets_ruv/difexp_none_ruv.tsv"
  ),
  list(
    id    = "2_3_all_datasets",
    label = "2nd Trim vs Term — All Datasets",
    combat = "output/phase2b_combat/phase2b_2_3_all_datasets/difexp_none_combat.tsv",
    ruv    = "output/phase2b_ruv/phase2b_2_3_all_datasets_ruv/difexp_none_ruv.tsv"
  )
)

out_dir <- "articles/imputation_article/misc"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (cmp in comparisons) {
  degs_combat <- load_degs(cmp$combat)
  degs_ruv    <- load_degs(cmp$ruv)

  cat(sprintf("\n%s:\n", cmp$label))
  cat(sprintf("  ComBat DEGs: %d\n", length(degs_combat)))
  cat(sprintf("  RUV DEGs: %d\n", length(degs_ruv)))
  cat(sprintf("  Overlap: %d\n", length(intersect(degs_combat, degs_ruv))))
  cat(sprintf("  ComBat only: %d\n", length(setdiff(degs_combat, degs_ruv))))
  cat(sprintf("  RUV only: %d\n", length(setdiff(degs_ruv, degs_combat))))

  out_file <- file.path(out_dir, sprintf("venn_ruv_vs_combat_%s.pdf", cmp$id))

  venn.plot <- venn.diagram(
    x = list(
      "ComBat" = degs_combat,
      "RUV (k=2)" = degs_ruv
    ),
    filename        = NULL,
    fill            = c("#4BACC6", "#9BBB59"),
    alpha           = 0.5,
    cat.cex         = 1.1,
    cat.fontface    = "bold",
    cex             = 1.4,
    fontface        = "bold",
    main            = cmp$label,
    main.cex        = 1.3,
    main.fontface   = "bold",
    sub             = sprintf("DEGs at FDR < %.2f (no logFC cutoff)", fdr_cutoff),
    sub.cex         = 0.9,
    margin          = 0.1
  )

  pdf(out_file, width = 6, height = 5)
  grid.draw(venn.plot)
  dev.off()
  cat(sprintf("  Wrote: %s\n", out_file))
}
