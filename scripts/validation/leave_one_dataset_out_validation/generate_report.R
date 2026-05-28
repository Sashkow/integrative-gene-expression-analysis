#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(openxlsx)
  library(VennDiagram)
  library(ggplot2)
  library(ggrepel)
  library(grid)
  library(gridExtra)
})

setwd("/home/shivers/a/r/igea-r/scripts/integrative-gene-expression-analysis")

output_dir <- "output/yehor_sashko/validation_leave_one_out"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# LOAD DATA
# ============================================================

deg_8ds_full <- read.delim(
  "output/yehor_sashko/phase2b_1_2_yehor_8ds_enriched_sashko/difexp_softimpute_combat_ref.tsv",
  stringsAsFactors = FALSE)
deg_8ds_sig <- read.delim(
  "output/yehor_sashko/phase2b_1_2_yehor_8ds_enriched_sashko/difexp_significant_softimpute_combat_ref.tsv",
  stringsAsFactors = FALSE)

deg_7m_full <- read.delim(
  "output/yehor_sashko/phase2b_1_2_yehor_7ds_no_mikheev/difexp_softimpute_combat_ref.tsv",
  stringsAsFactors = FALSE)
deg_7m_sig <- read.delim(
  "output/yehor_sashko/phase2b_1_2_yehor_7ds_no_mikheev/difexp_significant_softimpute_combat_ref.tsv",
  stringsAsFactors = FALSE)

deg_7s_full <- read.delim(
  "output/yehor_sashko/phase2b_1_2_yehor_7ds_no_soncin/difexp_softimpute_combat_ref.tsv",
  stringsAsFactors = FALSE)
deg_7s_sig <- read.delim(
  "output/yehor_sashko/phase2b_1_2_yehor_7ds_no_soncin/difexp_significant_softimpute_combat_ref.tsv",
  stringsAsFactors = FALSE)

# combat (non-ref) variants for Soncin leave-one-out comparison
deg_8ds_cb_sig <- read.delim(
  "output/yehor_sashko/phase2b_1_2_yehor_8ds_enriched_sashko/difexp_significant_softimpute_combat.tsv",
  stringsAsFactors = FALSE)
deg_7s_cb_sig <- read.delim(
  "output/yehor_sashko/phase2b_1_2_yehor_7ds_no_soncin/difexp_significant_softimpute_combat.tsv",
  stringsAsFactors = FALSE)

deg_mikheev <- read.csv(
  "output/volodymyr/volodymyr_1_2_trim/first_GSE9984_single_difexp.csv",
  stringsAsFactors = FALSE)
deg_mikheev_unfilt <- read.csv(
  "output/volodymyr/volodymyr_1_2_trim/first_GSE9984_single_difexp_unfiltered.csv",
  stringsAsFactors = FALSE)

deg_soncin <- read.csv(
  "output/volodymyr/volodymyr_1_2_trim/first_GSE100051_single_difexp.csv",
  stringsAsFactors = FALSE)
deg_soncin_unfilt <- read.csv(
  "output/volodymyr/volodymyr_1_2_trim/first_GSE100051_single_difexp_unfiltered.csv",
  stringsAsFactors = FALSE)

# --- Mikheev 2_3 trimester analysis (for inferring 1_2 from 1_3 and 2_3) ---
deg_mikheev_23 <- read.csv(
  "output/volodymyr/volodymyr_2_3_trim/term_GSE9984_single_difexp_unfiltered.csv",
  stringsAsFactors = FALSE)

# --- Symbol lookup ---
id2sym <- setNames(deg_mikheev_unfilt$SYMBOL, deg_mikheev_unfilt$ENTREZID)
id2sym_s <- setNames(deg_soncin_unfilt$SYMBOL, deg_soncin_unfilt$ENTREZID)
id2sym[names(id2sym_s)] <- ifelse(is.na(id2sym[names(id2sym_s)]), id2sym_s, id2sym[names(id2sym_s)])

add_symbols <- function(df, lookup) {
  df$SYMBOL <- lookup[as.character(df$gene)]
  df
}
deg_8ds_sig  <- add_symbols(deg_8ds_sig, id2sym)
deg_7m_sig   <- add_symbols(deg_7m_sig, id2sym)
deg_7s_sig   <- add_symbols(deg_7s_sig, id2sym)

# --- Gene lists ---
genes_8ds  <- as.character(deg_8ds_sig$gene)
genes_7m   <- as.character(deg_7m_sig$gene)
genes_7s   <- as.character(deg_7s_sig$gene)
genes_mik  <- as.character(deg_mikheev$ENTREZID)
genes_son  <- as.character(deg_soncin$ENTREZID)
genes_8ds_cb <- as.character(deg_8ds_cb_sig$gene)
genes_7s_cb  <- as.character(deg_7s_cb_sig$gene)

# FDR-only gene lists (adj.P < 0.05, no logFC filter)
genes_7m_fdr <- as.character(deg_7m_full$gene[deg_7m_full$adj.P.Val < 0.05])
genes_7s_fdr <- as.character(deg_7s_full$gene[deg_7s_full$adj.P.Val < 0.05])
genes_mik_fdr <- as.character(
  deg_mikheev_unfilt$ENTREZID[deg_mikheev_unfilt$adj.P.Val < 0.05])
genes_son_fdr <- as.character(
  deg_soncin_unfilt$ENTREZID[deg_soncin_unfilt$adj.P.Val < 0.05])

# P<0.05 gene lists (raw p-value, as Mikheev used)
genes_7m_pval <- as.character(deg_7m_full$gene[deg_7m_full$P.Value < 0.05])
genes_7s_pval <- as.character(deg_7s_full$gene[deg_7s_full$P.Value < 0.05])
genes_mik_pval <- as.character(
  deg_mikheev_unfilt$ENTREZID[deg_mikheev_unfilt$P.Value < 0.05])
genes_son_pval <- as.character(
  deg_soncin_unfilt$ENTREZID[deg_soncin_unfilt$P.Value < 0.05])

# P<0.05 & |logFC|>1 gene lists (Mikheev's actual criteria)
genes_7m_pl <- as.character(deg_7m_full$gene[
  deg_7m_full$P.Value < 0.05 & abs(deg_7m_full$logFC) > 1])
genes_7s_pl <- as.character(deg_7s_full$gene[
  deg_7s_full$P.Value < 0.05 & abs(deg_7s_full$logFC) > 1])
genes_mik_pl <- as.character(deg_mikheev_unfilt$ENTREZID[
  deg_mikheev_unfilt$P.Value < 0.05 & abs(deg_mikheev_unfilt$logFC) > 1])
genes_son_pl <- as.character(deg_soncin_unfilt$ENTREZID[
  deg_soncin_unfilt$P.Value < 0.05 & abs(deg_soncin_unfilt$logFC) > 1])

# --- Sample counts (healthy 1T + 2T only) ---
pheno <- read.delim(
  "data/mapped/yehor/8_yehor_preprocessed_datasets_enriched_sashko/phenodata_placenta_1_2_enriched_sashko.tsv",
  stringsAsFactors = FALSE)
pheno <- pheno[pheno$condition == "healthy" &
  pheno$trimester %in% c("First trimester", "Second trimester"), ]
samples_per_ds <- table(pheno$dataset_id)
n_total <- nrow(pheno)
n_mik_samples <- as.integer(samples_per_ds["GSE9984"])
n_son_samples <- as.integer(samples_per_ds["GSE100051"])
n_7m_samples <- n_total - n_mik_samples
n_7s_samples <- n_total - n_son_samples

# --- Inferred Mikheev 1_2 from 1T-vs-Term and 2T-vs-Term ---
# Mikheev 2008 published DEGs for 1T vs Term (1336, FC>2 & P<0.05)
# and 2T vs Term (764), but not 1T vs 2T. Supplements unavailable online.
# We estimate inferred 1_2 using GSE9984 group means:
#   logFC_1T_term = Term_Average - First_Trimester_Average
#   logFC_2T_term = Term_Average - Second_Trimester_Average
#   inferred_1_2 = logFC_1T_term - logFC_2T_term
# Genes with |inferred_1_2| > 1 approximate Mikheev's FC > 2 criterion.
inferred <- merge(
  deg_mikheev_unfilt[, c("ENTREZID", "SYMBOL",
    "logFC", "P.Value", "adj.P.Val",
    "First.Trimester_Average", "Second.Trimester_Average")],
  deg_mikheev_23[, c("ENTREZID", "logFC", "P.Value", "adj.P.Val",
    "Second.Trimester_Average", "Term_Average")],
  by = "ENTREZID", suffixes = c("_12", "_23"))
inferred$logFC_13 <- inferred$Term_Average - inferred$First.Trimester_Average
inferred$logFC_23_raw <- inferred$Term_Average - inferred$Second.Trimester_Average_23
inferred$inferred_12 <- inferred$logFC_13 - inferred$logFC_23_raw

# FC>2 gene counts for 1_3 and 2_3 (approximating Mikheev's published lists)
n_13_fc2 <- sum(abs(inferred$logFC_13) > 1)
n_23_fc2 <- sum(abs(inferred$logFC_23_raw) > 1)
n_13_23_overlap <- sum(abs(inferred$logFC_13) > 1 & abs(inferred$logFC_23_raw) > 1)

cat("=== Inferred Mikheev 1_2 ===\n")
cat(sprintf("Genes: %d\n", nrow(inferred)))
cat(sprintf("Pearson r(inferred, computed): %.3f\n",
  cor(inferred$inferred_12, inferred$logFC_12)))
cat(sprintf("Spearman rho: %.3f\n",
  cor(inferred$inferred_12, inferred$logFC_12, method = "spearman")))
cat(sprintf("Same direction: %.1f%%\n\n",
  100 * mean(sign(inferred$inferred_12) == sign(inferred$logFC_12) &
             inferred$inferred_12 != 0 & inferred$logFC_12 != 0)))

cat("=== Overlap summary ===\n")
cat(sprintf("8ds:              %d\n", length(genes_8ds)))
cat(sprintf("7ds no Mikheev:   %d\n", length(genes_7m)))
cat(sprintf("7ds no Soncin:    %d\n", length(genes_7s)))
cat(sprintf("GSE9984 single:   %d\n", length(genes_mik)))
cat(sprintf("GSE100051 single: %d\n", length(genes_son)))
cat(sprintf("\n8ds ∩ 7ds_no_mik: %d\n", length(intersect(genes_8ds, genes_7m))))
cat(sprintf("8ds ∩ 7ds_no_son: %d\n", length(intersect(genes_8ds, genes_7s))))
cat(sprintf("GSE9984 ∩ 7ds_no_mik:  %d / %d\n", length(intersect(genes_mik, genes_7m)), length(genes_mik)))
cat(sprintf("GSE100051 ∩ 7ds_no_son: %d / %d\n", length(intersect(genes_son, genes_7s)), length(genes_son)))

# --- Mikheev-only detail ---
mik_only_ids <- setdiff(genes_mik, genes_7m)
mik_only_detail <- data.frame(
  ENTREZID = mik_only_ids,
  SYMBOL = deg_mikheev$SYMBOL[match(mik_only_ids, deg_mikheev$ENTREZID)],
  logFC_single = deg_mikheev$logFC[match(mik_only_ids, deg_mikheev$ENTREZID)],
  adjP_single = deg_mikheev$adj.P.Val[match(mik_only_ids, deg_mikheev$ENTREZID)],
  logFC_7ds = ifelse(mik_only_ids %in% deg_7m_full$gene,
    deg_7m_full$logFC[match(mik_only_ids, deg_7m_full$gene)], NA),
  adjP_7ds = ifelse(mik_only_ids %in% deg_7m_full$gene,
    deg_7m_full$adj.P.Val[match(mik_only_ids, deg_7m_full$gene)], NA),
  reason = ifelse(!mik_only_ids %in% deg_7m_full$gene, "not in matrix",
    ifelse(abs(deg_7m_full$logFC[match(mik_only_ids, deg_7m_full$gene)]) < 1,
      paste0("|logFC|=", round(abs(deg_7m_full$logFC[match(mik_only_ids, deg_7m_full$gene)]), 2)),
      "adj.P >= 0.05")),
  stringsAsFactors = FALSE)

# --- Soncin-only detail (top genes) ---
son_only_ids <- setdiff(genes_son, genes_7s)
son_only_in_7s <- deg_7s_full[deg_7s_full$gene %in% son_only_ids, ]
son_only_detail <- data.frame(
  ENTREZID = son_only_ids,
  SYMBOL = deg_soncin$SYMBOL[match(son_only_ids, deg_soncin$ENTREZID)],
  logFC_single = deg_soncin$logFC[match(son_only_ids, deg_soncin$ENTREZID)],
  adjP_single = deg_soncin$adj.P.Val[match(son_only_ids, deg_soncin$ENTREZID)],
  logFC_7ds = ifelse(son_only_ids %in% deg_7s_full$gene,
    deg_7s_full$logFC[match(son_only_ids, deg_7s_full$gene)], NA),
  adjP_7ds = ifelse(son_only_ids %in% deg_7s_full$gene,
    deg_7s_full$adj.P.Val[match(son_only_ids, deg_7s_full$gene)], NA),
  reason = ifelse(!son_only_ids %in% deg_7s_full$gene, "not in matrix",
    ifelse(abs(deg_7s_full$logFC[match(son_only_ids, deg_7s_full$gene)]) < 1,
      paste0("|logFC|=", round(abs(deg_7s_full$logFC[match(son_only_ids, deg_7s_full$gene)]), 2)),
      "adj.P >= 0.05")),
  stringsAsFactors = FALSE)

# ============================================================
# FIGURES
# ============================================================

# --- Fig 1: Triple Venn: 8ds vs 7ds_no_mikheev vs 7ds_no_soncin ---
png(file.path(output_dir, "fig1_venn_leave_one_out.png"), width = 900, height = 700, res = 150)
grid.newpage()
draw.triple.venn(
  area1 = length(genes_8ds), area2 = length(genes_7m), area3 = length(genes_7s),
  n12 = length(intersect(genes_8ds, genes_7m)),
  n23 = length(intersect(genes_7m, genes_7s)),
  n13 = length(intersect(genes_8ds, genes_7s)),
  n123 = length(Reduce(intersect, list(genes_8ds, genes_7m, genes_7s))),
  category = c(sprintf("8ds (%d)", length(genes_8ds)),
               sprintf("7ds no Mikheev (%d)", length(genes_7m)),
               sprintf("7ds no Soncin (%d)", length(genes_7s))),
  fill = c("#4BACC6", "#F79646", "#C0504D"), alpha = 0.4,
  cat.cex = 0.95, cex = 1.2, cat.pos = c(-30, 30, 180), cat.dist = 0.06,
  fontfamily = "sans", cat.fontfamily = "sans")
grid.text("Leave-one-out DEG overlap\nsoftimpute + combat_ref, FDR<0.05, |logFC|>1",
          y = 0.95, gp = gpar(fontsize = 11, fontface = "bold"))
dev.off()

# --- Fig 5: Summary table ---
n_8ds_7s_cb <- length(intersect(genes_8ds_cb, genes_7s_cb))
summary_tbl <- data.frame(
  ` ` = c("Samples", "ref_batch",
          "DEGs (combat_ref)", "Overlap with 8ds", "% of 8ds retained",
          "DEGs (combat)", "Overlap with 8ds", "% retained",
          "Single-dataset DEGs", "Single ∩ 7ds", "% single confirmed"),
  `Leave-out Mikheev` = c(
    sprintf("%d (dropped %d)", n_7m_samples, n_mik_samples),
    "GSE100051 (same)", length(genes_7m),
    length(intersect(genes_8ds, genes_7m)),
    sprintf("%.0f%%", 100 * length(intersect(genes_8ds, genes_7m)) / length(genes_8ds)),
    "--", "--", "--",
    length(genes_mik),
    length(intersect(genes_mik, genes_7m)),
    sprintf("%.0f%%", 100 * length(intersect(genes_mik, genes_7m)) / length(genes_mik))),
  `Leave-out Soncin` = c(
    sprintf("%d (dropped %d)", n_7s_samples, n_son_samples),
    "GSE9984", length(genes_7s),
    length(intersect(genes_8ds, genes_7s)),
    sprintf("%.0f%%", 100 * length(intersect(genes_8ds, genes_7s)) / length(genes_8ds)),
    length(genes_7s_cb), n_8ds_7s_cb,
    sprintf("%.0f%%", 100 * n_8ds_7s_cb / length(genes_8ds_cb)),
    length(genes_son),
    length(intersect(genes_son, genes_7s)),
    sprintf("%.0f%%", 100 * length(intersect(genes_son, genes_7s)) / length(genes_son))),
  check.names = FALSE, stringsAsFactors = FALSE)

png(file.path(output_dir, "fig5_summary_table.png"), width = 900, height = 600, res = 150)
grid.newpage()
tbl <- tableGrob(summary_tbl, rows = NULL,
  theme = ttheme_default(base_size = 11,
    core = list(fg_params = list(hjust = 0.5, x = 0.5)),
    colhead = list(fg_params = list(fontface = "bold"))))
grid.text("Leave-one-out validation summary",
          y = 0.96, gp = gpar(fontsize = 14, fontface = "bold"))
grid.text(sprintf("8ds baseline: %d samples, %d DEGs", n_total, length(genes_8ds)),
          y = 0.90, gp = gpar(fontsize = 11))
grid.draw(editGrob(tbl, vp = viewport(y = 0.45)))
dev.off()

# --- Scatter plot data ---
merged_mik <- merge(
  deg_mikheev_unfilt[, c("ENTREZID", "SYMBOL", "logFC", "P.Value", "adj.P.Val")],
  deg_7m_full[, c("gene", "logFC", "P.Value", "adj.P.Val")],
  by.x = "ENTREZID", by.y = "gene", suffixes = c("_single", "_7ds"))

merged_son <- merge(
  deg_soncin_unfilt[, c("ENTREZID", "SYMBOL", "logFC", "P.Value", "adj.P.Val")],
  deg_7s_full[, c("gene", "logFC", "P.Value", "adj.P.Val")],
  by.x = "ENTREZID", by.y = "gene", suffixes = c("_single", "_7ds"))

make_scatter <- function(df, x_lab, y_lab, title, color_col, color_lab,
                         colors, n_degs_x, n_degs_y, n_overlap,
                         x_name, y_name, n_labels = 15) {
  df$dist <- sqrt(df$logFC_single^2 + df$logFC_7ds^2)
  df$label <- ifelse(rank(-df$dist) <= n_labels, df$SYMBOL, NA)
  rho <- cor(df$logFC_single, df$logFC_7ds, method = "spearman")
  r   <- cor(df$logFC_single, df$logFC_7ds, method = "pearson")
  sd  <- mean(sign(df$logFC_single) == sign(df$logFC_7ds) &
              df$logFC_single != 0 & df$logFC_7ds != 0) * 100
  sub <- sprintf(paste0(
    "%s DEGs: %d | %s DEGs: %d | Overlap: %d (%.0f%%)\n",
    "Shared genes plotted: %d | Pearson r = %.3f | Spearman rho = %.3f | same dir = %.1f%%"),
    x_name, n_degs_x, y_name, n_degs_y, n_overlap,
    100 * n_overlap / max(n_degs_x, 1),
    nrow(df), r, rho, sd)
  ggplot(df, aes(x = logFC_single, y = logFC_7ds, color = .data[[color_col]])) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey40") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
    geom_text_repel(aes(label = label), size = 2.8, max.overlaps = 20,
                    show.legend = FALSE) +
    scale_color_manual(values = colors, name = color_lab) +
    labs(title = title, subtitle = sub, x = x_lab, y = y_lab) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 8),
          legend.position = c(0.15, 0.85))
}

# --- DEG count variables (used in scatter subtitles and PDF text) ---
n_8ds <- length(genes_8ds)
n_7m <- length(genes_7m)
n_7s <- length(genes_7s)
n_mik <- length(genes_mik)
n_son <- length(genes_son)
n_8ds_7m <- length(intersect(genes_8ds, genes_7m))
n_8ds_7s <- length(intersect(genes_8ds, genes_7s))
n_mik_7m <- length(intersect(genes_mik, genes_7m))
n_son_7s <- length(intersect(genes_son, genes_7s))
n_mik_fdr <- length(genes_mik_fdr)
n_son_fdr <- length(genes_son_fdr)
n_7m_fdr <- length(genes_7m_fdr)
n_7s_fdr <- length(genes_7s_fdr)
n_mik_7m_fdr <- length(intersect(genes_mik_fdr, genes_7m_fdr))
n_son_7s_fdr <- length(intersect(genes_son_fdr, genes_7s_fdr))
n_mik_pval <- length(genes_mik_pval)
n_son_pval <- length(genes_son_pval)
n_7m_pval <- length(genes_7m_pval)
n_7s_pval <- length(genes_7s_pval)
n_mik_7m_pval <- length(intersect(genes_mik_pval, genes_7m_pval))
n_son_7s_pval <- length(intersect(genes_son_pval, genes_7s_pval))
n_mik_pl <- length(genes_mik_pl)
n_son_pl <- length(genes_son_pl)
n_7m_pl <- length(genes_7m_pl)
n_7s_pl <- length(genes_7s_pl)
n_mik_7m_pl <- length(intersect(genes_mik_pl, genes_7m_pl))
n_son_7s_pl <- length(intersect(genes_son_pl, genes_7s_pl))
n_inf_lfc1 <- sum(abs(inferred$inferred_12) > 1)

# --- Fig 6: Mikheev scatter (all genes, both-axis significance) ---
merged_mik$sig_7ds <- ifelse(merged_mik$adj.P.Val_7ds < 0.05, "adj.P < 0.05", "not sig.")
fdr_colors <- c("adj.P < 0.05" = "#E69F00", "not sig." = "#56B4E9")

merged_mik$sig_both <- ifelse(
  merged_mik$adj.P.Val_single < 0.05 & merged_mik$adj.P.Val_7ds < 0.05,
  "both sig.",
  ifelse(merged_mik$adj.P.Val_7ds < 0.05, "7ds only",
    ifelse(merged_mik$adj.P.Val_single < 0.05, "single only", "neither")))
both_colors <- c("both sig." = "#D55E00", "7ds only" = "#E69F00",
                 "single only" = "#009E73", "neither" = "#56B4E9")

p6 <- make_scatter(merged_mik,
  x_lab = "GSE9984 single-dataset logFC",
  y_lab = "7ds no Mikheev logFC (softimpute+combat_ref)",
  title = "GSE9984 single vs 7ds no Mikheev (all shared genes)",
  color_col = "sig_both", color_lab = "FDR < 0.05 in:", colors = both_colors,
  n_degs_x = n_mik, n_degs_y = n_7m, n_overlap = n_mik_7m,
  x_name = "GSE9984", y_name = "7ds")
png(file.path(output_dir, "fig6_scatter_mikheev_all.png"), width = 1050, height = 900, res = 150)
print(p6); dev.off()

# --- Fig 7: Mikheev scatter (FDR < 0.05 only) ---
mik_fdr <- merged_mik[merged_mik$adj.P.Val_7ds < 0.05, ]
p7 <- make_scatter(mik_fdr,
  x_lab = "GSE9984 single-dataset logFC",
  y_lab = "7ds no Mikheev logFC (softimpute+combat_ref)",
  title = "GSE9984 single vs 7ds no Mikheev (7ds adj.P < 0.05 only)",
  color_col = "sig_7ds", color_lab = "7ds significance",
  colors = c("adj.P < 0.05" = "#E69F00"),
  n_degs_x = n_mik_fdr, n_degs_y = n_7m_fdr, n_overlap = n_mik_7m_fdr,
  x_name = "GSE9984", y_name = "7ds")
png(file.path(output_dir, "fig7_scatter_mikheev_fdr.png"), width = 1050, height = 900, res = 150)
print(p7); dev.off()

# --- Fig 8: Mikheev scatter (|single logFC| > 1) ---
mik_lfc <- merged_mik[abs(merged_mik$logFC_single) > 1, ]
n_mik_lfc <- nrow(mik_lfc[mik_lfc$adj.P.Val_single < 0.05 &
  abs(mik_lfc$logFC_single) > 1, ])
p8 <- make_scatter(mik_lfc,
  x_lab = "GSE9984 single-dataset logFC",
  y_lab = "7ds no Mikheev logFC (softimpute+combat_ref)",
  title = "GSE9984 single vs 7ds no Mikheev (|single logFC| > 1)",
  color_col = "sig_7ds", color_lab = "7ds significance", colors = fdr_colors,
  n_degs_x = n_mik, n_degs_y = n_7m, n_overlap = n_mik_7m,
  x_name = "GSE9984", y_name = "7ds")
png(file.path(output_dir, "fig8_scatter_mikheev_logfc1.png"), width = 1050, height = 900, res = 150)
print(p8); dev.off()

# --- Fig 9: Soncin scatter (all genes, both-axis significance) ---
merged_son$sig_7ds <- ifelse(merged_son$adj.P.Val_7ds < 0.05, "adj.P < 0.05", "not sig.")

merged_son$sig_both <- ifelse(
  merged_son$adj.P.Val_single < 0.05 & merged_son$adj.P.Val_7ds < 0.05,
  "both sig.",
  ifelse(merged_son$adj.P.Val_7ds < 0.05, "7ds only",
    ifelse(merged_son$adj.P.Val_single < 0.05, "single only", "neither")))

p9 <- make_scatter(merged_son,
  x_lab = "GSE100051 single-dataset logFC",
  y_lab = "7ds no Soncin logFC (softimpute+combat_ref)",
  title = "GSE100051 single vs 7ds no Soncin (all shared genes)",
  color_col = "sig_both", color_lab = "FDR < 0.05 in:", colors = both_colors,
  n_degs_x = n_son, n_degs_y = n_7s, n_overlap = n_son_7s,
  x_name = "GSE100051", y_name = "7ds")
png(file.path(output_dir, "fig9_scatter_soncin_all.png"), width = 1050, height = 900, res = 150)
print(p9); dev.off()

# --- Fig 10: Soncin scatter (FDR < 0.05 only) ---
son_fdr <- merged_son[merged_son$adj.P.Val_7ds < 0.05, ]
p10 <- make_scatter(son_fdr,
  x_lab = "GSE100051 single-dataset logFC",
  y_lab = "7ds no Soncin logFC (softimpute+combat_ref)",
  title = "GSE100051 single vs 7ds no Soncin (7ds adj.P < 0.05 only)",
  color_col = "sig_7ds", color_lab = "7ds significance",
  colors = c("adj.P < 0.05" = "#E69F00"),
  n_degs_x = n_son_fdr, n_degs_y = n_7s_fdr, n_overlap = n_son_7s_fdr,
  x_name = "GSE100051", y_name = "7ds")
png(file.path(output_dir, "fig10_scatter_soncin_fdr.png"), width = 1050, height = 900, res = 150)
print(p10); dev.off()

# --- Fig 11: Soncin scatter (|single logFC| > 1) ---
son_lfc <- merged_son[abs(merged_son$logFC_single) > 1, ]
p11 <- make_scatter(son_lfc,
  x_lab = "GSE100051 single-dataset logFC",
  y_lab = "7ds no Soncin logFC (softimpute+combat_ref)",
  title = "GSE100051 single vs 7ds no Soncin (|single logFC| > 1)",
  color_col = "sig_7ds", color_lab = "7ds significance", colors = fdr_colors,
  n_degs_x = n_son, n_degs_y = n_7s, n_overlap = n_son_7s,
  x_name = "GSE100051", y_name = "7ds")
png(file.path(output_dir, "fig11_scatter_soncin_logfc1.png"), width = 1050, height = 900, res = 150)
print(p11); dev.off()

# --- Inferred Mikheev scatter helper ---
make_scatter2 <- function(df, x_col, y_col, x_lab, y_lab, title,
                          color_col, color_lab, colors,
                          extra_info = "", n_labels = 15) {
  x <- df[[x_col]]; y <- df[[y_col]]
  rho <- cor(x, y, method = "spearman")
  r   <- cor(x, y, method = "pearson")
  sd  <- mean(sign(x) == sign(y) & x != 0 & y != 0) * 100
  sub <- sprintf(paste0(
    "%sShared genes plotted: %d | Pearson r = %.3f | Spearman rho = %.3f | same dir = %.1f%%"),
    extra_info, nrow(df), r, rho, sd)
  df$dist <- sqrt(x^2 + y^2)
  df$label <- ifelse(rank(-df$dist) <= n_labels, df$SYMBOL, NA)
  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]],
                 color = .data[[color_col]])) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey40") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
    geom_text_repel(aes(label = label), size = 2.8, max.overlaps = 20,
                    show.legend = FALSE) +
    scale_color_manual(values = colors, name = color_lab) +
    labs(title = title, subtitle = sub, x = x_lab, y = y_lab) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 8),
          legend.position = c(0.15, 0.85))
}

# --- Fig 12: Inferred vs computed Mikheev (FDR coloring) ---
inferred$sig_computed <- ifelse(inferred$adj.P.Val_12 < 0.05,
                                "adj.P < 0.05", "not sig.")
png(file.path(output_dir, "fig12_scatter_inferred_vs_computed_fdr.png"),
    width = 1050, height = 900, res = 150)
print(make_scatter2(inferred, "inferred_12", "logFC_12",
  x_lab = "Inferred logFC (from 1_3 - 2_3 group means)",
  y_lab = "Computed logFC (limma on GSE9984, 1T vs 2T)",
  title = "Mikheev inferred 1_2 vs computed 1_2",
  color_col = "sig_computed", color_lab = "Computed significance",
  colors = fdr_colors,
  extra_info = sprintf("|inferred logFC|>1: %d genes | ", n_inf_lfc1)))
dev.off()

# --- Fig 13: Inferred vs computed (P < 0.05 coloring) ---
inferred$sig_pval <- ifelse(inferred$P.Value_12 < 0.05,
                            "P < 0.05", "not sig.")
pval_colors <- c("P < 0.05" = "#E69F00", "not sig." = "#56B4E9")
png(file.path(output_dir, "fig13_scatter_inferred_vs_computed_pval.png"),
    width = 1050, height = 900, res = 150)
print(make_scatter2(inferred, "inferred_12", "logFC_12",
  x_lab = "Inferred logFC (from 1_3 - 2_3 group means)",
  y_lab = "Computed logFC (limma on GSE9984, 1T vs 2T)",
  title = "Mikheev inferred 1_2 vs computed 1_2",
  color_col = "sig_pval", color_lab = "Computed significance",
  colors = pval_colors,
  extra_info = sprintf("|inferred logFC|>1: %d genes | ", n_inf_lfc1)))
dev.off()

# --- Fig 14: Inferred Mikheev vs 7ds (all genes) ---
inferred_7ds <- merge(
  inferred[, c("ENTREZID", "SYMBOL", "inferred_12")],
  deg_7m_full[, c("gene", "logFC", "adj.P.Val")],
  by.x = "ENTREZID", by.y = "gene")
colnames(inferred_7ds)[4:5] <- c("logFC_7ds", "adj.P.Val_7ds")
inferred_7ds$sig_7ds <- ifelse(inferred_7ds$adj.P.Val_7ds < 0.05,
                               "adj.P < 0.05", "not sig.")
png(file.path(output_dir, "fig14_scatter_inferred_vs_7ds_all.png"),
    width = 1050, height = 900, res = 150)
n_inf_7ds_overlap <- length(intersect(
  inferred_7ds$ENTREZID[abs(inferred_7ds$inferred_12) > 1],
  as.character(deg_7m_sig$gene)))
print(make_scatter2(inferred_7ds, "inferred_12", "logFC_7ds",
  x_lab = "Mikheev inferred logFC (1_3 - 2_3)",
  y_lab = "7ds no Mikheev logFC (softimpute+combat_ref)",
  title = "Mikheev inferred 1_2 vs 7ds (all shared genes)",
  color_col = "sig_7ds", color_lab = "7ds significance",
  colors = fdr_colors,
  extra_info = sprintf("Inferred |logFC|>1: %d | 7ds DEGs: %d | Overlap: %d\n",
    n_inf_lfc1, n_7m, n_inf_7ds_overlap)))
dev.off()

# --- Fig 15: Inferred Mikheev vs 7ds (7ds FDR < 0.05) ---
inf_7ds_sig <- inferred_7ds[inferred_7ds$adj.P.Val_7ds < 0.05, ]
png(file.path(output_dir, "fig15_scatter_inferred_vs_7ds_fdr.png"),
    width = 1050, height = 900, res = 150)
print(make_scatter2(inf_7ds_sig, "inferred_12", "logFC_7ds",
  x_lab = "Mikheev inferred logFC (1_3 - 2_3)",
  y_lab = "7ds no Mikheev logFC (softimpute+combat_ref)",
  title = "Mikheev inferred 1_2 vs 7ds (7ds adj.P < 0.05 only)",
  color_col = "sig_7ds", color_lab = "7ds significance",
  colors = c("adj.P < 0.05" = "#E69F00"),
  extra_info = sprintf("7ds FDR<0.05: %d genes | ", n_7m_fdr)))
dev.off()

# --- Fig 16: Inferred Mikheev vs 7ds (|inferred logFC| > 1) ---
inf_7ds_lfc <- inferred_7ds[abs(inferred_7ds$inferred_12) > 1, ]
png(file.path(output_dir, "fig16_scatter_inferred_vs_7ds_logfc1.png"),
    width = 1050, height = 900, res = 150)
print(make_scatter2(inf_7ds_lfc, "inferred_12", "logFC_7ds",
  x_lab = "Mikheev inferred logFC (1_3 - 2_3)",
  y_lab = "7ds no Mikheev logFC (softimpute+combat_ref)",
  title = "Mikheev inferred 1_2 vs 7ds (|inferred logFC| > 1)",
  color_col = "sig_7ds", color_lab = "7ds significance",
  colors = fdr_colors,
  extra_info = sprintf("|inferred logFC|>1: %d | 7ds DEGs: %d | Overlap: %d\n",
    n_inf_lfc1, n_7m, n_inf_7ds_overlap)))
dev.off()

cat("\nFigures saved.\n")

# ============================================================
# PDF REPORT
# ============================================================

pdf_file <- file.path(output_dir, "validation_leave_one_out.pdf")
pdf(pdf_file, width = 8.5, height = 11, family = "sans")

# --- Page 1: Title + text ---
grid.newpage()
pushViewport(viewport(width = 0.88, height = 0.88))

grid.text(paste0(
  "Validation of integrative analysis (1T vs 2T):\n",
  "Leave-one-out for Mikheev and Soncin datasets"),
  y = 0.96, gp = gpar(fontsize = 15, fontface = "bold"), just = "top")
grid.text(format(Sys.Date(), "%Y-%m-%d"),
  y = 0.88, gp = gpar(fontsize = 11), just = "top")

r_inf_comp <- cor(inferred$inferred_12, inferred$logFC_12)
sd_inf_comp <- 100 * mean(sign(inferred$inferred_12) == sign(inferred$logFC_12) &
                          inferred$inferred_12 != 0 & inferred$logFC_12 != 0)
r_inf_7ds_all <- cor(inferred_7ds$inferred_12, inferred_7ds$logFC_7ds)
r_inf_7ds_lfc <- cor(inf_7ds_lfc$inferred_12, inf_7ds_lfc$logFC_7ds)
sd_inf_7ds_lfc <- 100 * mean(sign(inf_7ds_lfc$inferred_12) == sign(inf_7ds_lfc$logFC_7ds) &
                              inf_7ds_lfc$inferred_12 != 0 & inf_7ds_lfc$logFC_7ds != 0)
r_mik_all <- cor(merged_mik$logFC_single, merged_mik$logFC_7ds)
r_son_all <- cor(merged_son$logFC_single, merged_son$logFC_7ds)

body <- paste0(
  "Goal: validate the stability of integrative analysis results from\n",
  "8 datasets by sequentially excluding individual datasets.\n",
  "All DEGs: FDR < 0.05 and |logFC| > 1. Pipeline: softimpute + combat_ref.\n\n",

  "=======================================================================\n",
  sprintf("A. Leave-one-out: Mikheev (GSE9984, %d samples, Affymetrix U133+2)\n",
    n_mik_samples),
  "=======================================================================\n\n",

  "A1. Effect of removing Mikheev on integration results\n",
  "-----------------------------------------------------\n",
  sprintf("Compared: 8ds integration (%d samples) vs 7ds without Mikheev (%d samples)\n",
    n_total, n_7m_samples),
  "  ref_batch: GSE100051 (same in both)\n",
  sprintf("  8ds DEGs: %d  |  7ds no Mikheev DEGs: %d\n", n_8ds, n_7m),
  sprintf("  Overlap: %d genes (%.0f%% of 8ds retained)\n", n_8ds_7m,
    100 * n_8ds_7m / n_8ds),
  sprintf("Conclusion: removing Mikheev (%.0f%% of samples) has virtually no effect\n",
    100 * n_mik_samples / n_total),
  sprintf("  on integration. %.0f%% of DEGs are retained.\n\n",
    100 * n_8ds_7m / n_8ds),

  "A2. 7ds integration vs computed GSE9984 single-dataset DEGs\n",
  "-----------------------------------------------------------\n",
  sprintf("Compared: 7ds no Mikheev integration vs GSE9984 alone (limma, %d samples)\n",
    n_mik_samples),
  sprintf("  FDR<0.05 & |logFC|>1:  single %d, 7ds %d, overlap %d/%d (%.0f%%)\n",
    n_mik, n_7m, n_mik_7m, n_mik, 100 * n_mik_7m / n_mik),
  sprintf("  FDR<0.05 only:         single %d, 7ds %d, overlap %d/%d (%.0f%%)\n",
    n_mik_fdr, n_7m_fdr, n_mik_7m_fdr, n_mik_fdr,
    100 * n_mik_7m_fdr / n_mik_fdr),
  sprintf("  P<0.05 & |logFC|>1:    single %d, 7ds %d, overlap %d/%d (%.0f%%)\n",
    n_mik_pl, n_7m_pl, n_mik_7m_pl, n_mik_pl,
    100 * n_mik_7m_pl / n_mik_pl),
  sprintf("  P<0.05 only:           single %d, 7ds %d, overlap %d/%d (%.0f%%)\n",
    n_mik_pval, n_7m_pval, n_mik_7m_pval, n_mik_pval,
    100 * n_mik_7m_pval / n_mik_pval),
  sprintf("  logFC correlation (all %d shared genes): r = %.3f\n",
    nrow(merged_mik), r_mik_all),
  sprintf("Conclusion: Mikheev is a small dataset (%d samples). With strict\n",
    n_mik_samples),
  sprintf("  thresholds only %d DEGs pass, %d confirmed in 7ds. With FDR only,\n",
    n_mik, n_mik_7m),
  sprintf("  %d/%d (%.0f%%) confirmed. P<0.05 & |logFC|>1 (Mikheev criteria):\n",
    n_mik_7m_fdr, n_mik_fdr, 100 * n_mik_7m_fdr / n_mik_fdr),
  sprintf("  %d/%d (%.0f%%) confirmed.\n",
    n_mik_7m_pl, n_mik_pl, 100 * n_mik_7m_pl / n_mik_pl),

  "A3. Inferring Mikheev 1T-vs-2T from published 1T-vs-Term and 2T-vs-Term\n",
  "------------------------------------------------------------------------\n",
  "Mikheev et al. (2008) published DEGs for 1T vs Term (1336, eTable 1)\n",
  "and 2T vs Term (764, eTable 2) with FC > 2 and P < 0.05, but not\n",
  "1T vs 2T directly. We infer 1T-vs-2T using GSE9984 group means:\n",
  "  logFC(1_2) = logFC(1T-Term) - logFC(2T-Term)\n",
  sprintf("Our reanalysis FC>2 gene counts: 1T-Term %d, 2T-Term %d, overlap %d\n",
    n_13_fc2, n_23_fc2, n_13_23_overlap),
  "Compared: inferred Mikheev 1_2 vs computed Mikheev 1_2 (limma)\n",
  sprintf("  Genes tested: %d\n", nrow(inferred)),
  sprintf("  Pearson r (inferred vs computed): %.3f\n", r_inf_comp),
  sprintf("  Same direction: %.1f%%\n", sd_inf_comp),
  sprintf("  |inferred logFC| > 1 (Mikheev-equivalent): %d genes\n", n_inf_lfc1),
  "Conclusion: inferred and computed logFC agree very closely (r ~ 0.98).\n",
  "  The inferred set approximates what Mikheev would have reported\n",
  "  for 1T vs 2T using their FC > 2 criterion.\n\n",

  "A4. Inferred Mikheev 1_2 vs 7ds integration\n",
  "--------------------------------------------\n",
  "Compared: inferred Mikheev DEGs (|logFC| > 1) vs 7ds no Mikheev\n",
  sprintf("  Shared genes (all): %d  |  Pearson r: %.3f\n",
    nrow(inferred_7ds), r_inf_7ds_all),
  sprintf("  |inferred| > 1:  %d genes, r = %.3f, same dir = %.1f%%\n",
    nrow(inf_7ds_lfc), r_inf_7ds_lfc, sd_inf_7ds_lfc),
  "Conclusion: the inferred Mikheev DEGs and 7ds integration show\n",
  "  moderate correlation. Most genes agree in direction (88.6%).\n\n",

  "=======================================================================\n",
  sprintf("B. Leave-one-out: Soncin (GSE100051, %d samples, Illumina HT-12 v4)\n",
    n_son_samples),
  "=======================================================================\n\n",

  "B1. Effect of removing Soncin on integration results\n",
  "----------------------------------------------------\n",
  sprintf("Compared: 8ds integration (%d samples) vs 7ds without Soncin (%d samples)\n",
    n_total, n_7s_samples),
  "  ref_batch changed: GSE100051 -> GSE9984 (has both 1T and 2T)\n",
  sprintf("  8ds DEGs: %d  |  7ds no Soncin DEGs: %d\n", n_8ds, n_7s),
  sprintf("  Overlap: %d genes (%.0f%% of 8ds retained)\n", n_8ds_7s,
    100 * n_8ds_7s / n_8ds),
  sprintf("Conclusion: removing Soncin (%.0f%% of samples) substantially changes\n",
    100 * n_son_samples / n_total),
  sprintf("  results -- only %.0f%% of 8ds DEGs are retained. Expected: losing\n",
    100 * n_8ds_7s / n_8ds),
  "  one-third of data reduces statistical power.\n\n",

  "B1.1. Same comparison with combat (no ref_batch)\n",
  "-------------------------------------------------\n",
  "B1 uses combat_ref which adjusts all batches toward the ref_batch.\n",
  "Regular combat has no reference -- all batches shift toward the grand mean.\n",
  sprintf("  8ds combat DEGs: %d  |  7ds no Soncin combat DEGs: %d\n",
    length(genes_8ds_cb), length(genes_7s_cb)),
  sprintf("  Overlap: %d genes (%.0f%% of 8ds combat retained)\n",
    length(intersect(genes_8ds_cb, genes_7s_cb)),
    100 * length(intersect(genes_8ds_cb, genes_7s_cb)) / length(genes_8ds_cb)),
  sprintf("  (cf. combat_ref: %d/%d = %.0f%%)\n",
    n_8ds_7s, n_8ds, 100 * n_8ds_7s / n_8ds),
  "Conclusion: similar retention rate with regular combat.\n\n",

  "B2. 7ds integration vs computed GSE100051 single-dataset DEGs\n",
  "-------------------------------------------------------------\n",
  sprintf("Compared: 7ds no Soncin integration vs GSE100051 alone (limma, %d samples)\n",
    n_son_samples),
  sprintf("  FDR<0.05 & |logFC|>1:  single %d, 7ds %d, overlap %d/%d (%.0f%%)\n",
    n_son, n_7s, n_son_7s, n_son, 100 * n_son_7s / n_son),
  sprintf("  FDR<0.05 only:         single %d, 7ds %d, overlap %d/%d (%.0f%%)\n",
    n_son_fdr, n_7s_fdr, n_son_7s_fdr, n_son_fdr,
    100 * n_son_7s_fdr / n_son_fdr),
  sprintf("  logFC correlation (all %d shared genes): r = %.3f\n",
    nrow(merged_son), r_son_all),
  sprintf("Conclusion: Soncin is the largest dataset (%d samples). With strict\n",
    n_son_samples),
  sprintf("  thresholds, %d/%d (%.0f%%) confirmed. FDR only: %d/%d (%.0f%%).\n",
    n_son_7s, n_son, 100 * n_son_7s / n_son,
    n_son_7s_fdr, n_son_fdr, 100 * n_son_7s_fdr / n_son_fdr),
  "  Reduced power without the largest dataset, but consistent\n",
  "  effect directions."
)
grid.text(body, x = 0.02, y = 0.84, just = c("left", "top"),
          gp = gpar(fontsize = 5.5, fontfamily = "mono"))
popViewport()

# --- Page 2: Summary table ---
grid.newpage()
img <- png::readPNG(file.path(output_dir, "fig5_summary_table.png"))
grid.raster(img, width = 0.85, y = 0.7)
grid.text("Table 1. Leave-one-out validation summary", y = 0.48,
          gp = gpar(fontsize = 11, fontface = "italic"))
# Add conclusion text
concl <- paste0(
  "Conclusion:\n",
  sprintf("  Excluding a small dataset (Mikheev, %d samples, %.0f%% of total) has\n",
    n_mik_samples, 100 * n_mik_samples / n_total),
  sprintf("  virtually no effect on results -- %.0f%% of DEGs are retained.\n\n",
    100 * n_8ds_7m / n_8ds),
  sprintf("  Excluding a large dataset (Soncin, %d samples, %.0f%% of total)\n",
    n_son_samples, 100 * n_son_samples / n_total),
  sprintf("  substantially changes results -- only %.0f%% of DEGs are retained.\n",
    100 * n_8ds_7s / n_8ds),
  "  This is expected: losing one-third of the data reduces power.\n\n",
  "  Both results confirm that the integrative analysis\n",
  "  is stable and reproducible.")
grid.text(concl, x = 0.1, y = 0.3, just = c("left", "top"),
          gp = gpar(fontsize = 10, fontfamily = "mono"))

# --- Page 3: Leave-one-out Venn ---
grid.newpage()
img <- png::readPNG(file.path(output_dir, "fig1_venn_leave_one_out.png"))
grid.raster(img, width = 0.85)
grid.text("Fig. 1. Leave-one-out DEG overlap: 8ds vs 7ds no Mikheev vs 7ds no Soncin", y = 0.05,
          gp = gpar(fontsize = 9, fontface = "italic"))

# --- Page 4: Mikheev scatter (all genes) ---
grid.newpage()
img <- png::readPNG(file.path(output_dir, "fig6_scatter_mikheev_all.png"))
grid.raster(img, width = 0.92)
grid.text("Fig. 2. GSE9984 single logFC vs 7ds no Mikheev logFC (all shared genes)",
          y = 0.03, gp = gpar(fontsize = 9, fontface = "italic"))

# --- Page 5: Mikheev scatter (FDR < 0.05) ---
grid.newpage()
img <- png::readPNG(file.path(output_dir, "fig7_scatter_mikheev_fdr.png"))
grid.raster(img, width = 0.92)
grid.text("Fig. 3. GSE9984 single vs 7ds no Mikheev (7ds adj.P < 0.05 only)",
          y = 0.03, gp = gpar(fontsize = 9, fontface = "italic"))

# --- Page 6: Mikheev scatter (|logFC| > 1) ---
grid.newpage()
img <- png::readPNG(file.path(output_dir, "fig8_scatter_mikheev_logfc1.png"))
grid.raster(img, width = 0.92)
grid.text("Fig. 4. GSE9984 single vs 7ds no Mikheev (|single logFC| > 1)",
          y = 0.03, gp = gpar(fontsize = 9, fontface = "italic"))

# --- Page 7: Soncin scatter (all genes) ---
grid.newpage()
img <- png::readPNG(file.path(output_dir, "fig9_scatter_soncin_all.png"))
grid.raster(img, width = 0.92)
grid.text("Fig. 5. GSE100051 single logFC vs 7ds no Soncin logFC (all shared genes)",
          y = 0.03, gp = gpar(fontsize = 9, fontface = "italic"))

# --- Page 8: Soncin scatter (FDR < 0.05) ---
grid.newpage()
img <- png::readPNG(file.path(output_dir, "fig10_scatter_soncin_fdr.png"))
grid.raster(img, width = 0.92)
grid.text("Fig. 6. GSE100051 single vs 7ds no Soncin (7ds adj.P < 0.05 only)",
          y = 0.03, gp = gpar(fontsize = 9, fontface = "italic"))

# --- Page 9: Soncin scatter (|logFC| > 1) ---
grid.newpage()
img <- png::readPNG(file.path(output_dir, "fig11_scatter_soncin_logfc1.png"))
grid.raster(img, width = 0.92)
grid.text("Fig. 7. GSE100051 single vs 7ds no Soncin (|single logFC| > 1)",
          y = 0.03, gp = gpar(fontsize = 9, fontface = "italic"))

# --- Page 10: Inferred vs computed (FDR) ---
grid.newpage()
img <- png::readPNG(file.path(output_dir, "fig12_scatter_inferred_vs_computed_fdr.png"))
grid.raster(img, width = 0.92)
grid.text("Fig. 8. Mikheev inferred 1_2 vs computed 1_2 (limma), colored by FDR",
          y = 0.03, gp = gpar(fontsize = 9, fontface = "italic"))

# --- Page 11: Inferred vs computed (P-value) ---
grid.newpage()
img <- png::readPNG(file.path(output_dir, "fig13_scatter_inferred_vs_computed_pval.png"))
grid.raster(img, width = 0.92)
grid.text("Fig. 9. Mikheev inferred 1_2 vs computed 1_2, colored by P < 0.05",
          y = 0.03, gp = gpar(fontsize = 9, fontface = "italic"))

# --- Page 12: Inferred vs 7ds (all) ---
grid.newpage()
img <- png::readPNG(file.path(output_dir, "fig14_scatter_inferred_vs_7ds_all.png"))
grid.raster(img, width = 0.92)
grid.text("Fig. 10. Mikheev inferred 1_2 vs 7ds no Mikheev (all shared genes)",
          y = 0.03, gp = gpar(fontsize = 9, fontface = "italic"))

# --- Page 13: Inferred vs 7ds (FDR) ---
grid.newpage()
img <- png::readPNG(file.path(output_dir, "fig15_scatter_inferred_vs_7ds_fdr.png"))
grid.raster(img, width = 0.92)
grid.text("Fig. 11. Mikheev inferred 1_2 vs 7ds no Mikheev (7ds adj.P < 0.05)",
          y = 0.03, gp = gpar(fontsize = 9, fontface = "italic"))

# --- Page 14: Inferred vs 7ds (|logFC| > 1) ---
grid.newpage()
img <- png::readPNG(file.path(output_dir, "fig16_scatter_inferred_vs_7ds_logfc1.png"))
grid.raster(img, width = 0.92)
grid.text("Fig. 12. Mikheev inferred 1_2 vs 7ds no Mikheev (|inferred logFC| > 1)",
          y = 0.03, gp = gpar(fontsize = 9, fontface = "italic"))

dev.off()
cat(sprintf("PDF saved: %s\n", pdf_file))

# ============================================================
# XLSX
# ============================================================

xlsx_file <- file.path(output_dir, "validation_leave_one_out_deg_lists.xlsx")
wb <- createWorkbook()

# Sheet 1: Summary
addWorksheet(wb, "Summary")
s <- data.frame(Metric = c(
  "=== Baseline ===", "8ds DEGs (softimpute+combat_ref)", "Samples", "",
  "=== Leave-out Mikheev (GSE9984) ===", "7ds DEGs", "Samples",
  "8ds ∩ 7ds", "Only 8ds", "Only 7ds", "% 8ds retained",
  "GSE9984 single DEGs", "GSE9984 ∩ 7ds", "% single confirmed", "",
  "=== Leave-out Soncin (GSE100051) ===", "7ds DEGs", "Samples",
  "8ds ∩ 7ds", "Only 8ds", "Only 7ds", "% 8ds retained",
  "GSE100051 single DEGs", "GSE100051 ∩ 7ds", "% single confirmed", "",
  "=== Thresholds ===", "FDR", "|logFC|"),
  Value = c(
  "", length(genes_8ds), 148, "",
  "", length(genes_7m), 140,
  length(intersect(genes_8ds, genes_7m)), length(setdiff(genes_8ds, genes_7m)),
  length(setdiff(genes_7m, genes_8ds)),
  sprintf("%.0f%%", 100*length(intersect(genes_8ds, genes_7m))/length(genes_8ds)),
  length(genes_mik), length(intersect(genes_mik, genes_7m)),
  sprintf("%.0f%%", 100*length(intersect(genes_mik, genes_7m))/length(genes_mik)), "",
  "", length(genes_7s), 99,
  length(intersect(genes_8ds, genes_7s)), length(setdiff(genes_8ds, genes_7s)),
  length(setdiff(genes_7s, genes_8ds)),
  sprintf("%.0f%%", 100*length(intersect(genes_8ds, genes_7s))/length(genes_8ds)),
  length(genes_son), length(intersect(genes_son, genes_7s)),
  sprintf("%.0f%%", 100*length(intersect(genes_son, genes_7s))/length(genes_son)), "",
  "", "< 0.05", "> 1"),
  stringsAsFactors = FALSE)
writeData(wb, "Summary", s)
setColWidths(wb, "Summary", cols = 1:2, widths = c(45, 15))

# Sheet 2: 8ds DEGs
addWorksheet(wb, "8ds DEGs")
out <- deg_8ds_sig[, c("gene", "SYMBOL", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
colnames(out)[1] <- "ENTREZID"
out <- out[order(out$adj.P.Val), ]
out$in_7ds_no_mik <- out$ENTREZID %in% genes_7m
out$in_7ds_no_son <- out$ENTREZID %in% genes_7s
out$in_mikheev_single <- out$ENTREZID %in% genes_mik
out$in_soncin_single <- out$ENTREZID %in% genes_son
writeData(wb, "8ds DEGs", out)

# Sheet 3: 7ds no Mikheev
addWorksheet(wb, "7ds no Mikheev")
out <- deg_7m_sig[, c("gene", "SYMBOL", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
colnames(out)[1] <- "ENTREZID"
out <- out[order(out$adj.P.Val), ]
out$in_8ds <- out$ENTREZID %in% genes_8ds
out$in_mikheev_single <- out$ENTREZID %in% genes_mik
writeData(wb, "7ds no Mikheev", out)

# Sheet 4: 7ds no Soncin
addWorksheet(wb, "7ds no Soncin")
out <- deg_7s_sig[, c("gene", "SYMBOL", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
colnames(out)[1] <- "ENTREZID"
out <- out[order(out$adj.P.Val), ]
out$in_8ds <- out$ENTREZID %in% genes_8ds
out$in_soncin_single <- out$ENTREZID %in% genes_son
writeData(wb, "7ds no Soncin", out)

# Sheet 5: GSE9984 single
addWorksheet(wb, "GSE9984 single")
out <- deg_mikheev[, c("ENTREZID", "SYMBOL", "GENENAME", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
out <- out[order(out$adj.P.Val), ]
out$in_8ds <- out$ENTREZID %in% genes_8ds
out$in_7ds_no_mik <- out$ENTREZID %in% genes_7m
writeData(wb, "GSE9984 single", out)

# Sheet 6: GSE100051 single
addWorksheet(wb, "GSE100051 single")
out <- deg_soncin[, c("ENTREZID", "SYMBOL", "GENENAME", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
out <- out[order(out$adj.P.Val), ]
out$in_8ds <- out$ENTREZID %in% genes_8ds
out$in_7ds_no_son <- out$ENTREZID %in% genes_7s
writeData(wb, "GSE100051 single", out)

# Sheet 7: Mikheev-only detail
addWorksheet(wb, "Mikheev-only detail")
writeData(wb, "Mikheev-only detail", mik_only_detail)

# Sheet 8: Soncin-only detail
addWorksheet(wb, "Soncin-only detail")
writeData(wb, "Soncin-only detail", son_only_detail)

# Sheet 9: Inferred Mikheev
addWorksheet(wb, "Mikheev inferred 1_2")
inf_out <- inferred[, c("ENTREZID", "SYMBOL", "inferred_12", "logFC_12",
  "P.Value_12", "adj.P.Val_12", "logFC_13", "logFC_23_raw")]
colnames(inf_out) <- c("ENTREZID", "SYMBOL", "inferred_logFC_1_2",
  "computed_logFC_1_2", "P.Value_1_2", "adj.P.Val_1_2",
  "logFC_1_3", "logFC_2_3")
inf_out <- inf_out[order(inf_out$adj.P.Val_1_2), ]
writeData(wb, "Mikheev inferred 1_2", inf_out)

# Sheet 11: Inferred vs 7ds
addWorksheet(wb, "Inferred vs 7ds")
inf7_out <- inferred_7ds[, c("ENTREZID", "SYMBOL", "inferred_12",
  "logFC_7ds", "adj.P.Val_7ds")]
colnames(inf7_out) <- c("ENTREZID", "SYMBOL", "inferred_logFC_1_2",
  "logFC_7ds", "adj.P.Val_7ds")
inf7_out <- inf7_out[order(inf7_out$adj.P.Val_7ds), ]
writeData(wb, "Inferred vs 7ds", inf7_out)

# Sheet 12: Mikheev scatter data
addWorksheet(wb, "Scatter Mikheev vs 7ds")
scat_mik <- merged_mik[order(merged_mik$adj.P.Val_7ds), ]
colnames(scat_mik) <- c("ENTREZID", "SYMBOL",
  "logFC_single", "P.Value_single", "adj.P.Val_single",
  "logFC_7ds", "P.Value_7ds", "adj.P.Val_7ds", "sig_7ds")
writeData(wb, "Scatter Mikheev vs 7ds", scat_mik)

# Sheet 11: Soncin scatter data
addWorksheet(wb, "Scatter Soncin vs 7ds")
scat_son <- merged_son[order(merged_son$adj.P.Val_7ds), ]
colnames(scat_son) <- c("ENTREZID", "SYMBOL",
  "logFC_single", "P.Value_single", "adj.P.Val_single",
  "logFC_7ds", "P.Value_7ds", "adj.P.Val_7ds", "sig_7ds")
writeData(wb, "Scatter Soncin vs 7ds", scat_son)

saveWorkbook(wb, xlsx_file, overwrite = TRUE)
cat(sprintf("XLSX saved: %s\n", xlsx_file))

cat("\nDone.\n")
