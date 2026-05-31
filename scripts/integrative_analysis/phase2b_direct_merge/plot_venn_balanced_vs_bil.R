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

contrasts <- list(
  list(
    id  = "1_2",
    label = "First vs. Second Trimester",
    bal = "output/phase2b_combat/phase2b_1_2_balanced/difexp_none_combat.tsv",
    bil = "output/phase2b_batch_in_limma/phase2b_1_2_all_datasets_batch_in_limma/difexp_none_batch_in_limma.tsv"
  ),
  list(
    id  = "2_3",
    label = "Second Trimester vs. Term",
    bal = "output/phase2b_combat/phase2b_2_3_balanced/difexp_none_combat.tsv",
    bil = "output/phase2b_batch_in_limma/phase2b_2_3_all_datasets_batch_in_limma/difexp_none_batch_in_limma.tsv"
  )
)

out_dir <- "articles/imputation_article/misc"

for (ctr in contrasts) {
  degs_bal <- load_degs(ctr$bal)
  degs_bil <- load_degs(ctr$bil)

  cat(sprintf("\n%s:\n", ctr$label))
  cat(sprintf("  Balanced DEGs: %d\n", length(degs_bal)))
  cat(sprintf("  Batch-in-limma DEGs: %d\n", length(degs_bil)))
  cat(sprintf("  Overlap: %d\n", length(intersect(degs_bal, degs_bil))))
  cat(sprintf("  Balanced only: %d\n", length(setdiff(degs_bal, degs_bil))))
  cat(sprintf("  Batch-in-limma only: %d\n", length(setdiff(degs_bil, degs_bal))))

  out_file <- file.path(out_dir, sprintf("venn_balanced_vs_bil_%s.pdf", ctr$id))

  venn.plot <- venn.diagram(
    x = list(
      "Balanced\n(ComBat)" = degs_bal,
      "Batch-in-limma\n(no ComBat)" = degs_bil
    ),
    filename        = NULL,
    fill            = c("#4BACC6", "#F79646"),
    alpha           = 0.5,
    cat.cex         = 1.1,
    cat.fontface    = "bold",
    cex             = 1.4,
    fontface        = "bold",
    main            = ctr$label,
    main.cex        = 1.4,
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
