#!/usr/bin/env Rscript
#
# Volcano plots for sex-stratified DE results (6ds run).
#
# Usage:
#   Rscript scripts/yehor_conference_2026/volcano_sex_stratified.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(org.Hs.eg.db)
  library(latex2exp)
  library(VennDiagram)
})

setwd("/home/shivers/a/r/igea-r/scripts/integrative-gene-expression-analysis")

data_dir <- paste0(
  "articles/yehor_conference_2026/data/",
  "phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko/",
  "sex_stratified"
)
output_dir <- file.path(data_dir, "plots")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

comparisons <- list(
  list(
    file = "difexp_1t_male_vs_female.tsv",
    title = "1T: Male vs Female",
    tag = "1t_male_vs_female"
  ),
  list(
    file = "difexp_2t_male_vs_female.tsv",
    title = "2T: Male vs Female",
    tag = "2t_male_vs_female"
  ),
  list(
    file = "difexp_males_1t_vs_2t.tsv",
    title = "Males: 2T vs 1T",
    tag = "males_1t_vs_2t"
  ),
  list(
    file = "difexp_females_1t_vs_2t.tsv",
    title = "Females: 2T vs 1T",
    tag = "females_1t_vs_2t"
  )
)

fdr_thresh <- 0.05
logfc_thresh <- 1.0
n_top_labels <- 20

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

make_volcano <- function(de, title, out_path) {
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

  de$pi_score <- abs(de$logFC) * de$neg_log10_fdr
  top <- de[de$status != "NS" & !is.na(de$SYMBOL), ]
  top <- top[order(-top$pi_score), ]
  top <- head(top, n_top_labels)

  cut_line <- -log10(fdr_thresh)

  n_up <- sum(de$status == "Up")
  n_down <- sum(de$status == "Down")
  subtitle <- sprintf(
    "%d up, %d down (FDR < %.2f, |logFC| > %.1f)",
    n_up, n_down, fdr_thresh, logfc_thresh
  )

  p <- ggplot(de, aes(x = logFC, y = neg_log10_fdr)) +
    geom_point(
      aes(color = status),
      alpha = 0.7, size = 1.2
    ) +
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
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "grey40")
    )

  if (nrow(top) > 0) {
    top_up <- top[top$logFC > 0, ]
    top_down <- top[top$logFC < 0, ]
    if (nrow(top_up) > 0) {
      p <- p + geom_label_repel(
        data = top_up,
        aes(label = SYMBOL),
        nudge_x = 1, direction = "y",
        force = 4, size = 3, max.overlaps = 20
      )
    }
    if (nrow(top_down) > 0) {
      p <- p + geom_label_repel(
        data = top_down,
        aes(label = SYMBOL),
        nudge_x = -1, direction = "y",
        force = 4, size = 3, max.overlaps = 20
      )
    }
  }

  png(out_path, width = 10, height = 8, units = "in",
      res = 150, type = "cairo")
  print(p)
  dev.off()
  cat(sprintf("  Saved: %s\n", out_path))
}

cat("=== Volcano plots for sex-stratified DE (6ds) ===\n\n")

for (comp in comparisons) {
  path <- file.path(data_dir, comp$file)
  cat(sprintf("Processing: %s\n", comp$title))
  de <- read.delim(path, stringsAsFactors = FALSE)
  de <- add_symbols(de)
  out_file <- file.path(output_dir, paste0("volcano_", comp$tag, ".png"))
  make_volcano(de, comp$title, out_file)
}

# --- Venn diagram: 1_2 vs 1_2_m vs 1_2_f ---
cat("\n=== Venn diagram: 1_2 vs 1_2_m vs 1_2_f ===\n")

base_dir_6ds <- paste0(
  "articles/yehor_conference_2026/data/",
  "phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko"
)

sig_all <- read.delim(
  file.path(base_dir_6ds, "difexp_significant_softimpute_combat_ref.tsv"),
  stringsAsFactors = FALSE
)
sig_males <- read.delim(
  file.path(data_dir, "difexp_significant_males_1t_vs_2t.tsv"),
  stringsAsFactors = FALSE
)
sig_females <- read.delim(
  file.path(data_dir, "difexp_significant_females_1t_vs_2t.tsv"),
  stringsAsFactors = FALSE
)

genes_all <- sig_all$gene
genes_m <- sig_males$gene
genes_f <- sig_females$gene

cat(sprintf("  1_2 (all):     %d DEGs\n", length(genes_all)))
cat(sprintf("  1_2_m (males): %d DEGs\n", length(genes_m)))
cat(sprintf("  1_2_f (females): %d DEGs\n", length(genes_f)))

male_unique <- setdiff(genes_m, union(genes_all, genes_f))
female_unique <- setdiff(genes_f, union(genes_all, genes_m))

de_m_full <- read.delim(
  file.path(data_dir, "difexp_males_1t_vs_2t.tsv"),
  stringsAsFactors = FALSE
)
de_f_full <- read.delim(
  file.path(data_dir, "difexp_females_1t_vs_2t.tsv"),
  stringsAsFactors = FALSE
)

mu_stats <- de_f_full[de_f_full$gene %in% male_unique, ]
mu_near_miss <- sum(
  mu_stats$adj.P.Val < 0.10 & abs(mu_stats$logFC) > 0.8,
  na.rm = TRUE
)
mu_concordant <- sum(
  sign(de_m_full$logFC[match(male_unique, de_m_full$gene)]) ==
  sign(de_f_full$logFC[match(male_unique, de_f_full$gene)]),
  na.rm = TRUE
)

fu_stats <- de_m_full[de_m_full$gene %in% female_unique, ]
fu_near_miss <- sum(
  fu_stats$adj.P.Val < 0.10 & abs(fu_stats$logFC) > 0.8,
  na.rm = TRUE
)
fu_concordant <- sum(
  sign(de_m_full$logFC[match(female_unique, de_m_full$gene)]) ==
  sign(de_f_full$logFC[match(female_unique, de_f_full$gene)]),
  na.rm = TRUE
)
fu_truly_unique <- length(female_unique) - fu_near_miss

venn_path <- file.path(output_dir, "venn_1_2_sex_stratified.png")

venn.plot <- venn.diagram(
  x = list(
    "1T vs 2T" = genes_all,
    "Males 1T vs 2T" = genes_m,
    "Females 1T vs 2T" = genes_f
  ),
  filename = NULL,
  fill = c("#8DA0CB", "#66C2A5", "#FC8D62"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.dist = c(0.08, 0.08, 0.06),
  margin = 0.15,
  main = "DEG overlap: 1T vs 2T comparisons",
  main.cex = 1.4
)
venn.plot <- venn.plot[!sapply(venn.plot, function(x) inherits(x, "rect"))]

png(venn_path, width = 14, height = 9, units = "in",
    res = 150, type = "cairo")
grid::grid.newpage()

grid::pushViewport(grid::viewport(
  x = 0.5, y = 0.5, width = 0.50, height = 0.90
))
grid::grid.draw(venn.plot)
grid::popViewport()

dev.off()
cat(sprintf("  Saved: %s\n", venn_path))

cat("\nDone.\n")
