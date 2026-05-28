#!/usr/bin/env Rscript
library(ggplot2)
library(ggrepel)

BASE <- "articles/yehor_conference_2026/data/phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko"
SEX_DIR <- file.path(BASE, "sex_stratified")
VENN_DIR <- file.path(SEX_DIR, "venn_m_vs_f")
OUTPUT <- "articles/yehor_conference_2026/data/scatterplot_logfc_m_vs_f_delta1.png"

de_m <- read.delim(file.path(SEX_DIR, "difexp_males_1t_vs_2t.tsv"))
de_f <- read.delim(file.path(SEX_DIR, "difexp_females_1t_vs_2t.tsv"))

sig_m_genes <- de_m$gene[de_m$adj.P.Val < 0.05]
sig_f_genes <- de_f$gene[de_f$adj.P.Val < 0.05]

shared <- intersect(de_m$gene, de_f$gene)

df <- data.frame(
  gene = shared,
  logfc_m = de_m$logFC[match(shared, de_m$gene)],
  logfc_f = de_f$logFC[match(shared, de_f$gene)],
  fdr_m = de_m$adj.P.Val[match(shared, de_m$gene)],
  fdr_f = de_f$adj.P.Val[match(shared, de_f$gene)],
  stringsAsFactors = FALSE
)

df$sig_m <- df$gene %in% sig_m_genes
df$sig_f <- df$gene %in% sig_f_genes

df$category <- "Not significant"
df$category[df$sig_m & df$sig_f] <- "Both (shared)"
df$category[df$sig_m & !df$sig_f] <- "Males only"
df$category[!df$sig_m & df$sig_f] <- "Females only"
df$category <- factor(df$category,
                      levels = c("Both (shared)", "Males only",
                                 "Females only", "Not significant"))

ann <- read.delim(file.path(VENN_DIR, "all_groups_annotated.tsv"))
gene2sym <- setNames(ann$SYMBOL, ann$ENTREZID)
df$symbol <- gene2sym[as.character(df$gene)]
df$symbol[is.na(df$symbol)] <- ""

df_sig <- df[df$category != "Not significant", ]

# Treat non-significant logFC as 0 for delta calculation
eff_m <- ifelse(df_sig$sig_m, df_sig$logfc_m, 0)
eff_f <- ifelse(df_sig$sig_f, df_sig$logfc_f, 0)
df_sig <- df_sig[abs(eff_m - eff_f) > 1, ]

df_sig$outside <- abs(df_sig$logfc_m) >= 1 | abs(df_sig$logfc_f) >= 1
df_sig$color_group <- paste0(
  df_sig$category,
  ifelse(df_sig$outside, "", " (small effect)"))

top <- df_sig[df_sig$outside, ]
top$dist <- sqrt(top$logfc_m^2 + top$logfc_f^2)
top_genes <- head(top[order(-top$dist), ], 15)
df_sig$label <- ifelse(df_sig$gene %in% top_genes$gene & df_sig$symbol != "",
                       df_sig$symbol, NA)

r_val <- cor(df$logfc_m, df$logfc_f)

cat_colors <- c(
  "Both (shared)" = "#D32F2F",
  "Males only" = "#1976D2",
  "Females only" = "#F57C00",
  "Both (shared) (small effect)" = "#F0AAAA",
  "Males only (small effect)" = "#A8CBE8",
  "Females only (small effect)" = "#FDCFA1")

df_sig$color_group <- factor(df_sig$color_group,
  levels = names(cat_colors))

p <- ggplot(df_sig, aes(x = logfc_m, y = logfc_f, color = color_group)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "darkgreen", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "darkgreen", alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey40", alpha = 0.5) +
  geom_hline(yintercept = c(-1, 1), linetype = "dotted",
             color = "steelblue", alpha = 0.4) +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted",
             color = "steelblue", alpha = 0.4) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_label_repel(aes(label = label), size = 2.5,
                   max.overlaps = 20, box.padding = 0.3,
                   show.legend = FALSE) +
  scale_color_manual(
    values = cat_colors,
    name = "FDR<0.05, |ΔlogFC|>1",
    breaks = c("Both (shared)", "Males only", "Females only",
               "Both (shared) (small effect)",
               "Males only (small effect)",
               "Females only (small effect)"),
    labels = c("Both (shared) |logFC|≥1",
               "Males only |logFC|≥1",
               "Females only |logFC|≥1",
               "Both (shared) |logFC|<1",
               "Males only |logFC|<1",
               "Females only |logFC|<1")) +
  labs(
    title = sprintf(
      "logFC: Males 1T→2T vs Females 1T→2T (r = %.3f, n = %d)",
      r_val, nrow(df_sig)),
    x = "Males logFC (1T vs 2T, 46+6 = 52 samples)",
    y = "Females logFC (1T vs 2T, 56+9 = 65 samples)"
  ) +
  scale_x_continuous(breaks = -6:6) +
  scale_y_continuous(breaks = -6:6) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = c(0.15, 0.85),
    legend.background = element_rect(fill = "white",
                                     color = "grey80"),
    plot.title = element_text(hjust = 0.5, face = "bold",
                              size = 12)
  ) +
  coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5))

png(OUTPUT, width = 8, height = 7, units = "in",
    res = 150, type = "cairo")
print(p)
dev.off()

cat("Saved:", OUTPUT, "\n")
cat(sprintf("Pearson r = %.3f, n = %d genes\n", r_val, nrow(df)))
cat(sprintf("Both: %d, Males only: %d, Females only: %d, Neither: %d\n",
            sum(df$category == "Both (shared)"),
            sum(df$category == "Males only"),
            sum(df$category == "Females only"),
            sum(df$category == "Not significant")))
