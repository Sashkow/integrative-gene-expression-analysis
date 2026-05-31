options(bitmapType = "cairo")
library(ggplot2)
library(ggrepel)
library(VennDiagram)
library(grid)

output_dir <- "output/article_validation/prior_study_concordance"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

base_dir <- "output/yehor_sashko/phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko"

# --- Load DEG tables ---

sig_softimpute <- read.delim(file.path(base_dir, "difexp_significant_softimpute_combat_ref.tsv"),
                             stringsAsFactors = FALSE)
sig_intersection <- read.delim(file.path(base_dir, "difexp_significant_none_combat_ref.tsv"),
                               stringsAsFactors = FALSE)

full_softimpute <- read.delim(file.path(base_dir, "difexp_softimpute_combat_ref.tsv"),
                              stringsAsFactors = FALSE)

lykhenko_all <- read.csv("output/dissertation/integrative_pipeline_disser/1_2/difexp/difexp_all.csv",
                         stringsAsFactors = FALSE)
lykhenko_sig <- read.csv("output/dissertation/integrative_pipeline_disser/1_2/difexp/difexp_filtered.csv",
                         stringsAsFactors = FALSE)

# --- Identify gained DEGs ---

gained_ids <- setdiff(sig_softimpute$gene, sig_intersection$gene)
lost_ids   <- setdiff(sig_intersection$gene, sig_softimpute$gene)

cat("softImpute DEGs:", nrow(sig_softimpute), "\n")
cat("Intersection DEGs:", nrow(sig_intersection), "\n")
cat("Gained DEGs:", length(gained_ids), "\n")
cat("Lost DEGs (sig in intersection only):", length(lost_ids), "\n")
cat("Shared DEGs:", length(intersect(sig_softimpute$gene, sig_intersection$gene)), "\n\n")

# --- Match gained DEGs to 2021 full limma ---

lykhenko_all$ENTREZID <- as.character(lykhenko_all$ENTREZID)
gained_in_2021 <- gained_ids[gained_ids %in% lykhenko_all$ENTREZID]
gained_absent  <- gained_ids[!gained_ids %in% lykhenko_all$ENTREZID]

cat("Gained DEGs found in 2021 table:", length(gained_in_2021),
    sprintf("(%.1f%%)\n", 100 * length(gained_in_2021) / length(gained_ids)))
cat("Gained DEGs NOT in 2021 table:", length(gained_absent),
    "(non-Affy platform genes)\n\n")

# --- Build concordance table ---

concordance <- data.frame(
  entrezid = gained_in_2021,
  stringsAsFactors = FALSE
)

idx_current <- match(concordance$entrezid, as.character(sig_softimpute$gene))
concordance$logFC_current <- sig_softimpute$logFC[idx_current]
concordance$fdr_current   <- sig_softimpute$adj.P.Val[idx_current]

idx_2021 <- match(concordance$entrezid, lykhenko_all$ENTREZID)
concordance$symbol     <- lykhenko_all$SYMBOL[idx_2021]
concordance$logFC_2021 <- lykhenko_all$logFC[idx_2021]
concordance$fdr_2021   <- lykhenko_all$adj.P.Val[idx_2021]

concordance$same_direction <- sign(concordance$logFC_current) == sign(concordance$logFC_2021)

concordance$category <- ifelse(
  concordance$same_direction & concordance$fdr_2021 < 0.05, "same_dir_significant",
  ifelse(concordance$same_direction & concordance$fdr_2021 < 0.20, "same_dir_near_sig",
  ifelse(concordance$same_direction, "same_dir_weak",
  "opposite_direction")))

cat("=== Concordance summary (gained DEGs in 2021 table) ===\n")
cat_table <- table(concordance$category)
for (nm in names(cat_table)) {
  cat(sprintf("  %-25s %3d (%.1f%%)\n", nm, cat_table[nm],
              100 * cat_table[nm] / nrow(concordance)))
}
cat(sprintf("\nTotal same direction: %d/%d (%.1f%%)\n",
            sum(concordance$same_direction), nrow(concordance),
            100 * sum(concordance$same_direction) / nrow(concordance)))
cat(sprintf("Already FDR<0.05 in 2021: %d/%d (%.1f%%)\n",
            sum(concordance$category == "same_dir_significant"), nrow(concordance),
            100 * sum(concordance$category == "same_dir_significant") / nrow(concordance)))

write.csv(concordance, file.path(output_dir, "gained_deg_concordance.csv"),
          row.names = FALSE)

# --- Figure 1: logFC scatter plot ---

concordance$fdr_2021_group <- ifelse(
  concordance$fdr_2021 < 0.05, "FDR < 0.05 in 2021",
  ifelse(concordance$fdr_2021 < 0.20, "FDR 0.05-0.20", "FDR > 0.20"))
concordance$fdr_2021_group <- factor(concordance$fdr_2021_group,
  levels = c("FDR < 0.05 in 2021", "FDR 0.05-0.20", "FDR > 0.20"))

r_val <- cor(concordance$logFC_current, concordance$logFC_2021, use = "complete.obs")

top_genes <- concordance[order(-abs(concordance$logFC_current)), ][1:10, ]

p_scatter <- ggplot(concordance, aes(x = logFC_2021, y = logFC_current, colour = fdr_2021_group)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_vline(xintercept = 0, colour = "grey80") +
  geom_text_repel(data = top_genes, aes(label = symbol), size = 3,
                  max.overlaps = 15, show.legend = FALSE) +
  scale_colour_manual(values = c("FDR < 0.05 in 2021" = "#2166AC",
                                  "FDR 0.05-0.20" = "#F4A582",
                                  "FDR > 0.20" = "#B2182B")) +
  annotate("text", x = min(concordance$logFC_2021, na.rm = TRUE),
           y = max(concordance$logFC_current, na.rm = TRUE),
           label = sprintf("r = %.3f\nn = %d", r_val, nrow(concordance)),
           hjust = 0, vjust = 1, size = 4) +
  labs(x = "logFC (Lykhenko 2021, 4 Affy datasets, 22 samples)",
       y = "logFC (this study, 7 datasets, 123 samples)",
       colour = "2021 significance",
       title = "Gained DEGs: logFC concordance with prior analysis") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "fig_logfc_scatter_gained.pdf"), p_scatter,
       width = 7, height = 6)
png(file.path(output_dir, "fig_logfc_scatter_gained.png"),
    width = 2100, height = 1800, res = 300, type = "cairo")
print(p_scatter)
invisible(dev.off())
cat("\nSaved logFC scatter plot\n")

# --- Figure 2: Volcano plot of 2021 data with gained DEGs highlighted ---

volcano_df <- data.frame(
  entrezid = as.character(lykhenko_all$ENTREZID),
  symbol   = lykhenko_all$SYMBOL,
  logFC    = lykhenko_all$logFC,
  fdr      = lykhenko_all$adj.P.Val,
  stringsAsFactors = FALSE
)
volcano_df$neg_log10_fdr <- -log10(pmax(volcano_df$fdr, 1e-50))
volcano_df$is_gained <- volcano_df$entrezid %in% gained_in_2021
volcano_df$group <- ifelse(volcano_df$is_gained, "Gained DEG", "Other")

volcano_df <- volcano_df[order(volcano_df$is_gained), ]

top_gained <- volcano_df[volcano_df$is_gained, ]
top_gained <- top_gained[order(-top_gained$neg_log10_fdr), ][1:10, ]

p_volcano <- ggplot(volcano_df, aes(x = logFC, y = neg_log10_fdr, colour = group)) +
  geom_point(data = volcano_df[!volcano_df$is_gained, ], alpha = 0.15, size = 0.8) +
  geom_point(data = volcano_df[volcano_df$is_gained, ], alpha = 0.8, size = 2) +
  geom_text_repel(data = top_gained, aes(label = symbol), size = 3,
                  max.overlaps = 15, show.legend = FALSE) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = c("Other" = "grey70", "Gained DEG" = "#D6604D")) +
  labs(x = expression(log[2]~"fold change (Lykhenko 2021)"),
       y = expression(-log[10]~"FDR"),
       colour = NULL,
       title = "Gained DEGs in the 2021 Lykhenko analysis") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "fig_volcano_2021_gained.pdf"), p_volcano,
       width = 7, height = 6)
png(file.path(output_dir, "fig_volcano_2021_gained.png"),
    width = 2100, height = 1800, res = 300, type = "cairo")
print(p_volcano)
invisible(dev.off())
cat("Saved volcano plot\n")

# --- Figure 3: Three-way Venn diagram ---

set_softimpute  <- as.character(sig_softimpute$gene)
set_intersection <- as.character(sig_intersection$gene)
set_lykhenko    <- as.character(lykhenko_sig$ENTREZID)

venn_list <- list(
  "softImpute (380)" = set_softimpute,
  "Intersection (182)" = set_intersection,
  "Lykhenko 2021 (327)" = set_lykhenko
)

png(file.path(output_dir, "fig_venn_three_way.png"), width = 2400, height = 2000, res = 300, type = "cairo")
venn.plot <- venn.diagram(
  x = venn_list,
  filename = NULL,
  fill = c("#2166AC", "#F4A582", "#4DAF4A"),
  alpha = 0.5,
  cex = 1.2,
  cat.cex = 0.9,
  cat.dist = 0.07,
  margin = 0.1,
  main = "DEG overlap across three analyses",
  main.cex = 1.1
)
grid.draw(venn.plot)
dev.off()

pdf(file.path(output_dir, "fig_venn_three_way.pdf"), width = 8, height = 7)
grid.draw(venn.plot)
dev.off()
cat("Saved Venn diagram\n")

# --- Summary output ---

summary_text <- sprintf("
=== Prior-study concordance summary ===

This study (7 datasets, 123 samples, softImpute + ComBat-ref):
  Total DEGs: %d
  Intersection-only DEGs: %d
  Gained DEGs: %d (lost: %d)

Gained DEGs matched to Lykhenko 2021 (4 Affy datasets, 22 samples):
  Found in 2021 table: %d/%d (%.1f%%)
  Not found (non-Affy platform genes): %d

Direction concordance (of %d matched):
  Same direction, FDR<0.05 in 2021: %d (%.1f%%)
  Same direction, FDR 0.05-0.20:    %d (%.1f%%)
  Same direction, FDR>0.20:         %d (%.1f%%)
  Opposite direction:                %d (%.1f%%)

  Total same direction: %d/%d (%.1f%%)
  Pearson r of logFC:   %.3f

Three-way Venn set sizes:
  softImpute DEGs: %d
  Intersection DEGs: %d
  Lykhenko 2021 DEGs: %d
  softImpute ∩ Lykhenko 2021: %d
  Intersection ∩ Lykhenko 2021: %d
  All three: %d
",
nrow(sig_softimpute), nrow(sig_intersection), length(gained_ids), length(lost_ids),
length(gained_in_2021), length(gained_ids),
100 * length(gained_in_2021) / length(gained_ids),
length(gained_absent),
nrow(concordance),
sum(concordance$category == "same_dir_significant"),
100 * sum(concordance$category == "same_dir_significant") / nrow(concordance),
sum(concordance$category == "same_dir_near_sig"),
100 * sum(concordance$category == "same_dir_near_sig") / nrow(concordance),
sum(concordance$category == "same_dir_weak"),
100 * sum(concordance$category == "same_dir_weak") / nrow(concordance),
sum(concordance$category == "opposite_direction"),
100 * sum(concordance$category == "opposite_direction") / nrow(concordance),
sum(concordance$same_direction), nrow(concordance),
100 * sum(concordance$same_direction) / nrow(concordance),
r_val,
length(set_softimpute), length(set_intersection), length(set_lykhenko),
length(intersect(set_softimpute, set_lykhenko)),
length(intersect(set_intersection, set_lykhenko)),
length(Reduce(intersect, venn_list))
)

cat(summary_text)
writeLines(summary_text, file.path(output_dir, "concordance_summary.txt"))
cat("\nAll outputs saved to:", output_dir, "\n")
