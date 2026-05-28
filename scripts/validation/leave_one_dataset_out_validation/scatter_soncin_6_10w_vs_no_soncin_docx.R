#!/usr/bin/env Rscript
# Side-by-side scatters for Soncin leave-one-out (using Soncin 6-10w)
# A: GSE100051 (6-10w) single vs 7ds_no_soncin
# B: 8ds vs 7ds_no_soncin

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(officer)
})

setwd("/home/shivers/a/r/igea-r/scripts/integrative-gene-expression-analysis")

output_dir <- "output/yehor_sashko/validation_leave_one_out"

# ============================================================
# Data
# ============================================================

deg_son_unfilt <- read.csv(
  "output/yehor_sashko/soncin_6_10w_single/GSE100051_6_10w_single_difexp_unfiltered.csv",
  stringsAsFactors = FALSE)
deg_8ds_full <- read.delim(
  "output/yehor_sashko/phase2b_1_2_yehor_8ds_soncin_6_10w/difexp_softimpute_combat_ref.tsv",
  stringsAsFactors = FALSE)
deg_7s_full <- read.delim(
  "output/yehor_sashko/phase2b_1_2_yehor_7ds_no_soncin/difexp_softimpute_combat_ref.tsv",
  stringsAsFactors = FALSE)

# symbol lookup
id2sym <- setNames(deg_son_unfilt$SYMBOL, deg_son_unfilt$ENTREZID)

# sample counts
pheno <- read.delim(
  "data/mapped/yehor/8_yehor_preprocessed_datasets_enriched_sashko/phenodata_placenta_1_2_enriched_sashko.tsv",
  stringsAsFactors = FALSE)
pheno <- pheno[pheno$condition == "healthy" &
  pheno$trimester %in% c("First trimester", "Second trimester"), ]
son_pheno <- pheno[pheno$dataset_id == "GSE100051", ]
n_son_6_10w <- sum(
  (son_pheno$trimester == "First trimester" &
   son_pheno$fetus_week >= 6 & son_pheno$fetus_week <= 10) |
  son_pheno$trimester == "Second trimester")
n_8ds <- nrow(pheno) - nrow(son_pheno) + n_son_6_10w
n_7s  <- nrow(pheno) - nrow(son_pheno)

# --- Plot A: GSE100051 (6-10w) single vs 7ds_no_soncin ---
merged_a <- merge(
  deg_son_unfilt[, c("ENTREZID", "SYMBOL", "logFC", "P.Value", "adj.P.Val")],
  deg_7s_full[, c("gene", "logFC", "P.Value", "adj.P.Val")],
  by.x = "ENTREZID", by.y = "gene", suffixes = c("_single", "_7s"))

n_all_a <- nrow(merged_a)
df_a <- merged_a[merged_a$adj.P.Val_single < 0.05 |
                 merged_a$adj.P.Val_7s < 0.05, ]

df_a$sig_cat <- ifelse(
  df_a$adj.P.Val_single < 0.05 & df_a$adj.P.Val_7s < 0.05,
  "both FDR < 0.05",
  ifelse(df_a$adj.P.Val_7s < 0.05, "7ds only", "GSE100051 only"))

r_a   <- cor(df_a$logFC_single, df_a$logFC_7s, method = "pearson")
rho_a <- cor(df_a$logFC_single, df_a$logFC_7s, method = "spearman")
sd_a  <- mean(sign(df_a$logFC_single) == sign(df_a$logFC_7s) &
              df_a$logFC_single != 0 & df_a$logFC_7s != 0) * 100
n_both_a <- sum(df_a$sig_cat == "both FDR < 0.05")
n_7s_a   <- sum(df_a$sig_cat == "7ds only")
n_sing_a <- sum(df_a$sig_cat == "GSE100051 only")

df_a$dist <- sqrt(df_a$logFC_single^2 + df_a$logFC_7s^2)
df_a$label <- ifelse(rank(-df_a$dist) <= 15, df_a$SYMBOL, NA)

# --- Plot B: 8ds vs 7ds_no_soncin ---
merged_b <- merge(
  deg_8ds_full[, c("gene", "logFC", "P.Value", "adj.P.Val")],
  deg_7s_full[, c("gene", "logFC", "P.Value", "adj.P.Val")],
  by = "gene", suffixes = c("_8ds", "_7s"))
merged_b$SYMBOL <- id2sym[as.character(merged_b$gene)]

n_all_b <- nrow(merged_b)
df_b <- merged_b[merged_b$adj.P.Val_8ds < 0.05 |
                 merged_b$adj.P.Val_7s < 0.05, ]

df_b$sig_cat <- ifelse(
  df_b$adj.P.Val_8ds < 0.05 & df_b$adj.P.Val_7s < 0.05,
  "both FDR < 0.05",
  ifelse(df_b$adj.P.Val_7s < 0.05, "7ds only", "8ds only"))

r_b   <- cor(df_b$logFC_8ds, df_b$logFC_7s, method = "pearson")
rho_b <- cor(df_b$logFC_8ds, df_b$logFC_7s, method = "spearman")
sd_b  <- mean(sign(df_b$logFC_8ds) == sign(df_b$logFC_7s) &
              df_b$logFC_8ds != 0 & df_b$logFC_7s != 0) * 100
n_both_b <- sum(df_b$sig_cat == "both FDR < 0.05")
n_7s_b   <- sum(df_b$sig_cat == "7ds only")
n_8ds_b  <- sum(df_b$sig_cat == "8ds only")

df_b$dist <- sqrt(df_b$logFC_8ds^2 + df_b$logFC_7s^2)
df_b$label <- ifelse(rank(-df_b$dist) <= 15, df_b$SYMBOL, NA)

# ============================================================
# Plots
# ============================================================

sig_colors_a <- c("both FDR < 0.05" = "#E69F00",
                  "7ds only" = "#56B4E9",
                  "GSE100051 only" = "#009E73")

pa <- ggplot(df_a, aes(x = logFC_single, y = logFC_7s, color = sig_cat)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey40") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20,
                  show.legend = FALSE) +
  scale_color_manual(values = sig_colors_a, name = NULL) +
  labs(title = "A. GSE100051 single vs 7ds no Soncin",
       x = sprintf("GSE100051 logFC (%d samples)", n_son_6_10w),
       y = sprintf("7ds no Soncin logFC (%d samples)", n_7s)) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 13),
        legend.position = "inside",
        legend.position.inside = c(0.24, 0.88),
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.35, "cm"),
        legend.background = element_rect(fill = alpha("white", 0.7),
                                         color = NA),
        plot.margin = margin(5, 10, 5, 5))

sig_colors_b <- c("both FDR < 0.05" = "#E69F00",
                  "7ds only" = "#009E73",
                  "8ds only" = "#56B4E9")

pb <- ggplot(df_b, aes(x = logFC_8ds, y = logFC_7s, color = sig_cat)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey40") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20,
                  show.legend = FALSE) +
  scale_color_manual(values = sig_colors_b, name = NULL) +
  labs(title = "B. 8ds vs 7ds no Soncin",
       x = sprintf("8ds logFC (%d samples)", n_8ds),
       y = sprintf("7ds no Soncin logFC (%d samples)", n_7s)) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 13),
        legend.position = "inside",
        legend.position.inside = c(0.24, 0.88),
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.35, "cm"),
        legend.background = element_rect(fill = alpha("white", 0.7),
                                         color = NA),
        plot.margin = margin(5, 5, 5, 10))

combined <- pa + pb

img_path <- file.path(output_dir,
  "scatter_soncin_6_10w_vs_no_soncin_side_by_side.png")
png(img_path, width = 2400, height = 1100, res = 150)
print(combined)
dev.off()

# ============================================================
# DOCX
# ============================================================

deg_son_sig <- read.csv(
  "output/yehor_sashko/soncin_6_10w_single/GSE100051_6_10w_single_difexp.csv",
  stringsAsFactors = FALSE)
deg_8ds_sig <- read.delim(
  "output/yehor_sashko/phase2b_1_2_yehor_8ds_soncin_6_10w/difexp_significant_softimpute_combat_ref.tsv",
  stringsAsFactors = FALSE)
deg_7s_sig <- read.delim(
  "output/yehor_sashko/phase2b_1_2_yehor_7ds_no_soncin/difexp_significant_softimpute_combat_ref.tsv",
  stringsAsFactors = FALSE)

n_deg_son <- nrow(deg_son_sig)
n_deg_8ds <- nrow(deg_8ds_sig)
n_deg_7s  <- nrow(deg_7s_sig)
n_son_7s_overlap <- length(intersect(
  as.character(deg_son_sig$ENTREZID), as.character(deg_7s_sig$gene)))
n_8ds_7s_overlap <- length(intersect(
  as.character(deg_8ds_sig$gene), as.character(deg_7s_sig$gene)))

caption_a <- sprintf(paste0(
  "A. GSE100051 single-dataset, 1T restricted to weeks 6-10 (limma, %d samples, %d DEGs) vs ",
  "7ds no Soncin (%d samples, %d DEGs). ",
  "DEG overlap: %d/%d (%.0f%%). ",
  "Genes plotted: %d / %d (at least one FDR < 0.05). ",
  "Both sig.: %d, 7ds only: %d, GSE100051 only: %d. ",
  "Pearson r = %.3f, Spearman rho = %.3f, same direction = %.1f%%."),
  n_son_6_10w, n_deg_son, n_7s, n_deg_7s,
  n_son_7s_overlap, n_deg_son,
  100 * n_son_7s_overlap / max(n_deg_son, 1),
  nrow(df_a), n_all_a, n_both_a, n_7s_a, n_sing_a, r_a, rho_a, sd_a)

caption_b <- sprintf(paste0(
  "B. 8ds (%d samples, %d DEGs) vs 7ds no Soncin (%d samples, %d DEGs). ",
  "DEG overlap: %d (%.0f%% of 8ds retained). ",
  "Genes plotted: %d / %d (at least one FDR < 0.05). ",
  "Both sig.: %d, 7ds only: %d, 8ds only: %d. ",
  "Pearson r = %.3f, Spearman rho = %.3f, same direction = %.1f%%."),
  n_8ds, n_deg_8ds, n_7s, n_deg_7s,
  n_8ds_7s_overlap, 100 * n_8ds_7s_overlap / max(n_deg_8ds, 1),
  nrow(df_b), n_all_b, n_both_b, n_7s_b, n_8ds_b, r_b, rho_b, sd_b)

caption <- paste0(
  "Figure. Leave-one-out validation: effect of removing Soncin (GSE100051). ",
  "Soncin first-trimester samples restricted to gestational weeks 6-10 ",
  "(weeks 4, 5, 11, 12 excluded). ",
  "Only genes with FDR < 0.05 in at least one comparison are shown. ",
  "Dotted line = identity (y = x). DEG thresholds: FDR < 0.05 and |logFC| > 1. ",
  "Pipeline: softimpute + combat_ref.\n\n",
  caption_a, "\n\n", caption_b)

docx_path <- file.path(output_dir, "scatter_soncin_6_10w_vs_no_soncin.docx")
doc <- read_docx()
doc <- body_add_img(doc, img_path, width = 7.5, height = 3.4)
doc <- body_add_par(doc, caption, style = "Normal")
print(doc, target = docx_path)
cat(sprintf("DOCX saved: %s\n", docx_path))
