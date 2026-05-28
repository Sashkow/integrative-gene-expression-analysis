#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
})

setwd("/home/shivers/a/r/igea-r/scripts/integrative-gene-expression-analysis")

datasets <- c("GSE100051", "GSE93520", "GSE28551")

for (ds in datasets) {
  cat(sprintf("\n=== GO enrichment: CVS vs %s ===\n", ds))

  out_dir <- sprintf("output/cvs_vs_termination_compare/cvs_vs_%s", ds)
  de_file <- file.path(out_dir, "difexp_significant_none_ruv.tsv")

  if (!file.exists(de_file)) {
    cat(sprintf("  Skipping — %s not found\n", de_file))
    next
  }

  de <- read.delim(de_file, stringsAsFactors = FALSE)
  cat(sprintf("  Significant DEGs: %d (up=%d, down=%d)\n",
              nrow(de), sum(de$logFC > 0), sum(de$logFC < 0)))

  if (nrow(de) < 5) {
    cat("  Too few DEGs for enrichment, skipping\n")
    write.csv(data.frame(Note = "Too few DEGs for enrichment"),
              file.path(out_dir, "go_bp_up_in_cvs.csv"), row.names = FALSE)
    write.csv(data.frame(Note = "Too few DEGs for enrichment"),
              file.path(out_dir, "go_bp_down_in_cvs.csv"), row.names = FALSE)
    next
  }

  all_genes <- read.delim(file.path(out_dir, "difexp_none_ruv.tsv"),
                          stringsAsFactors = FALSE)
  universe <- rownames(all_genes)
  if (is.null(universe) || length(universe) == 0)
    universe <- all_genes[, 1]

  for (direction in c("up", "down")) {
    if (direction == "up") {
      genes <- rownames(de[de$logFC > 0, ])
      if (is.null(genes) || length(genes) == 0)
        genes <- de[de$logFC > 0, 1]
      label <- "Up in CVS"
    } else {
      genes <- rownames(de[de$logFC < 0, ])
      if (is.null(genes) || length(genes) == 0)
        genes <- de[de$logFC < 0, 1]
      label <- "Down in CVS"
    }

    if (length(genes) < 3) {
      cat(sprintf("  %s: only %d genes, skipping\n", label, length(genes)))
      next
    }

    ego <- enrichGO(gene      = genes,
                    universe  = universe,
                    OrgDb     = org.Hs.eg.db,
                    ont       = "BP",
                    keyType   = "ENTREZID",
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
                    minGSSize = 10,
                    maxGSSize = 500,
                    readable  = TRUE)

    csv_file <- file.path(out_dir, sprintf("go_bp_%s_in_cvs.csv", direction))
    if (!is.null(ego) && nrow(ego@result[ego@result$p.adjust < 0.05, ]) > 0) {
      res <- ego@result[ego@result$p.adjust < 0.05, ]
      write.csv(res, csv_file, row.names = FALSE)
      cat(sprintf("  %s: %d significant GO terms -> %s\n",
                  label, nrow(res), csv_file))

      png_file <- file.path(out_dir, sprintf("go_bp_%s_in_cvs_dotplot.png", direction))
      tryCatch({
        png(png_file, width = 10, height = 8, units = "in", res = 200)
        p <- dotplot(ego, showCategory = min(20, nrow(res)),
                     title = sprintf("GO BP %s — CVS vs %s", label, ds))
        print(p)
        dev.off()
        cat(sprintf("  Dotplot: %s\n", png_file))
      }, error = function(e) {
        try(dev.off(), silent = TRUE)
        cat(sprintf("  Dotplot failed: %s\n", e$message))
      })
    } else {
      write.csv(data.frame(Note = "No significant GO terms at FDR < 0.05"),
                csv_file, row.names = FALSE)
      cat(sprintf("  %s: no significant GO terms\n", label))
    }
  }
}
cat("\nGO enrichment complete.\n")
