#!/usr/bin/env Rscript
library(ggplot2)
library(ggrepel)
library(openxlsx)

BASE <- "articles/yehor_conference_2026/data/phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko"
SEX_DIR <- file.path(BASE, "sex_stratified")
VENN_DIR <- file.path(SEX_DIR, "venn_m_vs_f")
PRATER_PATH <- "articles/yehor_conference_2026/references/Prater_2021_RNA-Seq_first-second_trimester_transition_supp_tables.xlsx"
OUTPUT <- "articles/yehor_conference_2026/data/scatterplot_logfc_ours_vs_prater.png"

ours <- read.delim(file.path(BASE, "difexp_softimpute_combat_ref.tsv"))
sig_ours <- read.delim(file.path(BASE, "difexp_significant_softimpute_combat_ref.tsv"))

prater_wb <- loadWorkbook(PRATER_PATH)
prater <- read.xlsx(prater_wb, sheet = "T1 DEGs_results_table_l2fc1")
prater$entrez <- as.character(prater$entrez)
prater <- prater[!is.na(prater$entrez) & prater$entrez != "NA", ]

ours_genes <- setNames(ours$logFC, ours$gene)
ours_fdr <- setNames(ours$adj.P.Val, ours$gene)
ours_sig_set <- names(ours_fdr)[!is.na(ours_fdr) & ours_fdr < 0.05]

prater_logfc <- setNames(as.numeric(prater$log2FoldChange), prater$entrez)
prater_padj <- setNames(as.numeric(prater$padj), prater$entrez)
prater_sig <- prater$entrez[!is.na(prater_padj) & prater_padj < 0.05]

shared <- intersect(names(ours_genes), names(prater_logfc))

df <- data.frame(
  gene = shared,
  ours_logfc = ours_genes[shared],
  prater_logfc = prater_logfc[shared],
  stringsAsFactors = FALSE
)
df$ours_sig <- df$gene %in% ours_sig_set
df$prater_sig <- df$gene %in% prater_sig

df$category <- "Not significant"
df$category[df$ours_sig & df$prater_sig] <- "Both significant"
df$category[df$ours_sig & !df$prater_sig] <- "Our 1_2 only"
df$category[!df$ours_sig & df$prater_sig] <- "Prater only"
df$category <- factor(df$category,
                      levels = c("Both significant", "Our 1_2 only",
                                 "Prater only", "Not significant"))

ann <- read.delim(file.path(VENN_DIR, "all_groups_annotated.tsv"))
gene2sym <- setNames(ann$SYMBOL, ann$ENTREZID)
df$symbol <- gene2sym[df$gene]
df$symbol[is.na(df$symbol)] <- ""

both_sig <- df[df$category == "Both significant", ]
both_sig$dist <- sqrt(both_sig$ours_logfc^2 + both_sig$prater_logfc^2)
top_genes <- head(both_sig[order(-both_sig$dist), ], 15)
df$label <- ifelse(df$gene %in% top_genes$gene & df$symbol != "", df$symbol, NA)

r_val <- cor(df$ours_logfc, df$prater_logfc)

cat_colors <- c("Both significant" = "#D32F2F", "Our 1_2 only" = "#1976D2",
                "Prater only" = "#F57C00", "Not significant" = "#BDBDBD")

p <- ggplot(df, aes(x = ours_logfc, y = prater_logfc, color = category)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgreen", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgreen", alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40", alpha = 0.5) +
  geom_hline(yintercept = c(-1, 1), linetype = "dotted", color = "steelblue", alpha = 0.4) +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "steelblue", alpha = 0.4) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_label_repel(aes(label = label), size = 2.5, max.overlaps = 20,
                   box.padding = 0.3, show.legend = FALSE) +
  scale_color_manual(values = cat_colors, name = "Significance (FDR<0.05)") +
  labs(
    title = sprintf("logFC comparison: Our 6ds vs Prater 2021 (r = %.3f, n = %d)",
                    r_val, nrow(df)),
    x = "Our 6ds logFC (1T vs 2T, combined, 4–12 wk vs 14–19 wk, n=117)",
    y = "Prater 2021 log2FC (7–8 wk vs 13–14 wk, n=14)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = c(0.15, 0.85),
    legend.background = element_rect(fill = "white", color = "grey80"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12)
  ) +
  scale_x_continuous(breaks = -4:4) +
  scale_y_continuous(breaks = -5:5) +
  coord_cartesian(xlim = c(-4, 4), ylim = c(-5, 5))

png(OUTPUT, width = 8, height = 7, units = "in", res = 150, type = "cairo")
print(p)
dev.off()

cat("Saved:", OUTPUT, "\n")
cat(sprintf("Pearson r = %.3f, n = %d shared genes\n", r_val, nrow(df)))
cat(sprintf("Both sig: %d, Ours only: %d, Prater only: %d, Neither: %d\n",
            sum(df$category == "Both significant"),
            sum(df$category == "Our 1_2 only"),
            sum(df$category == "Prater only"),
            sum(df$category == "Not significant")))
