#!/usr/bin/env Rscript
# Generate 3 full-slide volcano plots: combined, males, females 1T vs 2T
# with top 10 gene labels each.

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(org.Hs.eg.db)
  library(latex2exp)
})

setwd("/home/shivers/a/r/igea-r/scripts/integrative-gene-expression-analysis")

base_dir <- paste0(
  "articles/yehor_conference_2026/data/",
  "phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko"
)
sex_dir <- file.path(base_dir, "sex_stratified")
output_dir <- file.path(sex_dir, "plots")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

fdr_thresh <- 0.05
logfc_thresh <- 1.0
n_top_labels <- 10

comparisons <- list(
  list(
    file = file.path(base_dir, "difexp_softimpute_combat_ref.tsv"),
    title = "Combined: 2T vs 1T (447 DEGs)",
    tag = "combined_1t_vs_2t"
  ),
  list(
    file = file.path(sex_dir, "difexp_males_1t_vs_2t.tsv"),
    title = "Males: 2T vs 1T (337 DEGs)",
    tag = "males_1t_vs_2t_top10",
    force = c("LEP")
  ),
  list(
    file = file.path(sex_dir, "difexp_females_1t_vs_2t.tsv"),
    title = "Females: 2T vs 1T (504 DEGs)",
    tag = "females_1t_vs_2t_top10"
  )
)

add_symbols <- function(de) {
  entrez_ids <- as.character(de$gene)
  ann <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = entrez_ids,
    columns = "SYMBOL",
    keytype = "ENTREZID"
  )
  ann <- ann[!duplicated(ann$ENTREZID), ]
  de$SYMBOL <- ann$SYMBOL[match(entrez_ids, ann$ENTREZID)]
  de
}

make_volcano <- function(de, title, out_path, force_symbols = NULL) {
  de$neg_log10_fdr <- -log10(de$adj.P.Val)
  de$neg_log10_fdr[is.infinite(de$neg_log10_fdr)] <-
    max(de$neg_log10_fdr[is.finite(de$neg_log10_fdr)]) * 1.05

  de$status <- ifelse(
    de$adj.P.Val < fdr_thresh & de$logFC > logfc_thresh, "Up",
    ifelse(
      de$adj.P.Val < fdr_thresh & de$logFC < -logfc_thresh, "Down",
      "NS"
    )
  )

  top <- de[de$status != "NS" & !is.na(de$SYMBOL), ]
  top <- top[order(-abs(top$logFC)), ]
  top <- head(top, n_top_labels)

  if (!is.null(force_symbols)) {
    forced <- de[de$SYMBOL %in% force_symbols & !is.na(de$SYMBOL), ]
    forced <- forced[!forced$SYMBOL %in% top$SYMBOL, ]
    top <- rbind(top, forced)
  }

  cut_line <- -log10(fdr_thresh)

  n_up <- sum(de$status == "Up")
  n_down <- sum(de$status == "Down")
  subtitle <- sprintf(
    "%d up, %d down (FDR < %.2f, |logFC| > %.1f)",
    n_up, n_down, fdr_thresh, logfc_thresh
  )

  p <- ggplot(de, aes(x = logFC, y = neg_log10_fdr)) +
    geom_point(aes(color = status), alpha = 0.7, size = 1.8) +
    scale_color_manual(
      values = c("Down" = "blue", "NS" = "grey60", "Up" = "red"),
      labels = c("Downregulated", "Not significant", "Upregulated")
    ) +
    geom_vline(
      xintercept = c(-logfc_thresh, logfc_thresh),
      linetype = "dashed", color = "red", linewidth = 0.5
    ) +
    geom_hline(
      yintercept = cut_line,
      linetype = "dashed", color = "red", linewidth = 0.5
    ) +
    labs(
      x = TeX("$\\log_{2}$ Fold Change"),
      y = TeX("$-\\log_{10}$ (Adjusted P-Value)"),
      title = title,
      subtitle = subtitle
    ) +
    theme_minimal(base_size = 16) +
    theme(
      legend.position = "top",
      legend.text = element_text(size = 14),
      legend.title = element_blank(),
      axis.title = element_text(size = 18, face = "bold"),
      plot.title = element_text(size = 22, face = "bold"),
      plot.subtitle = element_text(size = 14, color = "grey40")
    )

  if (nrow(top) > 0) {
    top_up <- top[top$logFC > 0, ]
    top_down <- top[top$logFC < 0, ]
    if (nrow(top_up) > 0) {
      p <- p + geom_label_repel(
        data = top_up, aes(label = SYMBOL),
        nudge_x = 1, direction = "y",
        force = 4, size = 4.5, max.overlaps = 20
      )
    }
    if (nrow(top_down) > 0) {
      p <- p + geom_label_repel(
        data = top_down, aes(label = SYMBOL),
        nudge_x = -1, direction = "y",
        force = 4, size = 4.5, max.overlaps = 20
      )
    }
  }

  png(out_path, width = 10, height = 7.5, units = "in",
      res = 150, type = "cairo")
  print(p)
  dev.off()
  cat(sprintf("  Saved: %s\n", out_path))
}

cat("=== Generating trio volcano plots (top 10 labels) ===\n\n")

for (comp in comparisons) {
  cat(sprintf("Processing: %s\n", comp$title))
  de <- read.delim(comp$file, stringsAsFactors = FALSE)
  de <- add_symbols(de)
  out_file <- file.path(output_dir, paste0("volcano_", comp$tag, ".png"))
  make_volcano(de, comp$title, out_file, force_symbols = comp$force)
}

cat("\nDone.\n")
