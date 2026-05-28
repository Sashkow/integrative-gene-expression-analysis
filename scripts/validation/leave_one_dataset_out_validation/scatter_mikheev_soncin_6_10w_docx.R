#!/usr/bin/env Rscript
# Mikheev leave-one-out side-by-side, using Soncin 6-10w filtered dataset
# A: GSE9984 single vs 7ds_no_mikheev_soncin_6_10w
# B: 7ds_no_mikheev_soncin_6_10w vs 8ds_soncin_6_10w

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

deg_mikheev_unfilt <- read.csv(
  "output/volodymyr/volodymyr_1_2_trim/first_GSE9984_single_difexp_unfiltered.csv",
  stringsAsFactors = FALSE)
deg_7m_full <- read.delim(
  "output/yehor_sashko/phase2b_1_2_yehor_7ds_no_mikheev_soncin_6_10w/difexp_softimpute_combat_ref.tsv",
  stringsAsFactors = FALSE)
deg_8ds_full <- read.delim(
  "output/yehor_sashko/phase2b_1_2_yehor_8ds_soncin_6_10w/difexp_softimpute_combat_ref.tsv",
  stringsAsFactors = FALSE)

# sample counts
pheno <- read.delim(
  "data/mapped/yehor/8_yehor_preprocessed_datasets_enriched_sashko/phenodata_placenta_1_2_enriched_sashko.tsv",
  stringsAsFactors = FALSE)
pheno <- pheno[pheno$condition == "healthy" &
  pheno$trimester %in% c("First trimester", "Second trimester"), ]
son_pheno <- pheno[pheno$dataset_id == "GSE100051", ]
n_son_kept <- sum(
  (son_pheno$trimester == "First trimester" &
   son_pheno$fetus_week >= 6 & son_pheno$fetus_week <= 10) |
  son_pheno$trimester == "Second trimester")
n_mik <- sum(pheno$dataset_id == "GSE9984")
n_8ds_filt <- nrow(pheno) - nrow(son_pheno) + n_son_kept
n_7m_filt <- n_8ds_filt - n_mik

# --- Plot A: GSE9984 single vs 7ds_no_mikheev_soncin_6_10w ---
merged_a <- merge(
  deg_mikheev_unfilt[, c("ENTREZID", "SYMBOL", "logFC", "P.Value", "adj.P.Val")],
  deg_7m_full[, c("gene", "logFC", "P.Value", "adj.P.Val")],
  by.x = "ENTREZID", by.y = "gene", suffixes = c("_single", "_7ds"))

n_all_a <- nrow(merged_a)
df_a <- merged_a[merged_a$adj.P.Val_single < 0.05 |
                 merged_a$adj.P.Val_7ds < 0.05, ]

df_a$sig_cat <- ifelse(
  df_a$adj.P.Val_single < 0.05 & df_a$adj.P.Val_7ds < 0.05,
  "both FDR < 0.05",
  ifelse(df_a$adj.P.Val_7ds < 0.05, "7ds only", "GSE9984 only"))

r_a   <- cor(df_a$logFC_single, df_a$logFC_7ds, method = "pearson")
rho_a <- cor(df_a$logFC_single, df_a$logFC_7ds, method = "spearman")
sd_a  <- mean(sign(df_a$logFC_single) == sign(df_a$logFC_7ds) &
              df_a$logFC_single != 0 & df_a$logFC_7ds != 0) * 100
n_both_a <- sum(df_a$sig_cat == "both FDR < 0.05")
n_7ds_a  <- sum(df_a$sig_cat == "7ds only")
n_sing_a <- sum(df_a$sig_cat == "GSE9984 only")

df_a$dist <- sqrt(df_a$logFC_single^2 + df_a$logFC_7ds^2)
df_a$label <- ifelse(rank(-df_a$dist) <= 15, df_a$SYMBOL, NA)

# --- Plot B: 7ds_no_mikheev_soncin_6_10w vs 8ds_soncin_6_10w ---
merged_b <- merge(
  deg_7m_full[, c("gene", "logFC", "P.Value", "adj.P.Val")],
  deg_8ds_full[, c("gene", "logFC", "P.Value", "adj.P.Val")],
  by = "gene", suffixes = c("_7ds", "_8ds"))
merged_b$SYMBOL <- deg_mikheev_unfilt$SYMBOL[
  match(merged_b$gene, deg_mikheev_unfilt$ENTREZID)]

n_all_b <- nrow(merged_b)
df_b <- merged_b[merged_b$adj.P.Val_7ds < 0.05 |
                 merged_b$adj.P.Val_8ds < 0.05, ]

df_b$sig_cat <- ifelse(
  df_b$adj.P.Val_7ds < 0.05 & df_b$adj.P.Val_8ds < 0.05,
  "both FDR < 0.05",
  ifelse(df_b$adj.P.Val_8ds < 0.05, "8ds only", "7ds only"))

r_b   <- cor(df_b$logFC_7ds, df_b$logFC_8ds, method = "pearson")
rho_b <- cor(df_b$logFC_7ds, df_b$logFC_8ds, method = "spearman")
sd_b  <- mean(sign(df_b$logFC_7ds) == sign(df_b$logFC_8ds) &
              df_b$logFC_7ds != 0 & df_b$logFC_8ds != 0) * 100
n_both_b <- sum(df_b$sig_cat == "both FDR < 0.05")
n_8ds_b  <- sum(df_b$sig_cat == "8ds only")
n_7ds_b  <- sum(df_b$sig_cat == "7ds only")

df_b$dist <- sqrt(df_b$logFC_7ds^2 + df_b$logFC_8ds^2)
df_b$label <- ifelse(rank(-df_b$dist) <= 15, df_b$SYMBOL, NA)

# ============================================================
# Plots
# ============================================================

sig_colors_a <- c("both FDR < 0.05" = "#E69F00",
                  "7ds only" = "#56B4E9",
                  "GSE9984 only" = "#009E73")

pa <- ggplot(df_a, aes(x = logFC_single, y = logFC_7ds, color = sig_cat)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey40") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20,
                  show.legend = FALSE) +
  scale_color_manual(values = sig_colors_a, name = NULL) +
  labs(title = "A. GSE9984 vs 7ds no Mikheev (Soncin 6-10w)",
       x = sprintf("GSE9984 logFC (%d samples)", n_mik),
       y = sprintf("7ds no Mikheev logFC (%d samples)", n_7m_filt)) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 12),
        legend.position = "inside",
        legend.position.inside = c(0.22, 0.88),
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.35, "cm"),
        legend.background = element_rect(fill = alpha("white", 0.7),
                                         color = NA),
        plot.margin = margin(5, 10, 5, 5))

sig_colors_b <- c("both FDR < 0.05" = "#E69F00",
                  "8ds only" = "#009E73",
                  "7ds only" = "#56B4E9")

pb <- ggplot(df_b, aes(x = logFC_7ds, y = logFC_8ds, color = sig_cat)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey40") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20,
                  show.legend = FALSE) +
  scale_color_manual(values = sig_colors_b, name = NULL) +
  labs(title = "B. 7ds no Mikheev vs 8ds (Soncin 6-10w)",
       x = sprintf("7ds no Mikheev logFC (%d samples)", n_7m_filt),
       y = sprintf("8ds logFC (%d samples)", n_8ds_filt)) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 12),
        legend.position = "inside",
        legend.position.inside = c(0.22, 0.88),
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.35, "cm"),
        legend.background = element_rect(fill = alpha("white", 0.7),
                                         color = NA),
        plot.margin = margin(5, 5, 5, 10))

combined <- pa + pb

img_path <- file.path(output_dir,
  "scatter_mikheev_7ds_8ds_soncin_6_10w_side_by_side.png")
png(img_path, width = 2400, height = 1100, res = 150)
print(combined)
dev.off()

# ============================================================
# DOCX
# ============================================================

deg_7m_sig <- read.delim(
  "output/yehor_sashko/phase2b_1_2_yehor_7ds_no_mikheev_soncin_6_10w/difexp_significant_softimpute_combat_ref.tsv",
  stringsAsFactors = FALSE)
deg_8ds_sig <- read.delim(
  "output/yehor_sashko/phase2b_1_2_yehor_8ds_soncin_6_10w/difexp_significant_softimpute_combat_ref.tsv",
  stringsAsFactors = FALSE)
deg_mik_sig <- read.csv(
  "output/volodymyr/volodymyr_1_2_trim/first_GSE9984_single_difexp.csv",
  stringsAsFactors = FALSE)

n_deg_7m <- nrow(deg_7m_sig)
n_deg_8ds <- nrow(deg_8ds_sig)
n_deg_mik <- nrow(deg_mik_sig)
n_mik_7m_overlap <- length(intersect(
  as.character(deg_mik_sig$ENTREZID), as.character(deg_7m_sig$gene)))
n_7m_8ds_overlap <- length(intersect(
  as.character(deg_7m_sig$gene), as.character(deg_8ds_sig$gene)))

caption_a <- sprintf(paste0(
  "Scatterplot A compares the Mikheev dataset (GSE9984) analyzed alone ",
  "(%d samples, %d DEGs) against the integrated dataset of the ",
  "remaining 7 datasets (%d samples, %d DEGs). ",
  "Of the %d Mikheev DEGs, %d (%.0f%%) were also identified in the 7-dataset integration. ",
  "Among %d genes plotted (out of %d shared, at least one with FDR < 0.05), ",
  "%d were significant in both, %d in the 7-dataset integration only, ",
  "and %d in Mikheev only. ",
  "Log-fold-change agreement: Pearson r = %.3f, Spearman rho = %.3f, ",
  "same direction in %.1f%% of genes."),
  n_mik, n_deg_mik, n_7m_filt, n_deg_7m,
  n_deg_mik, n_mik_7m_overlap,
  100 * n_mik_7m_overlap / max(n_deg_mik, 1),
  nrow(df_a), n_all_a, n_both_a, n_7ds_a, n_sing_a, r_a, rho_a, sd_a)

caption_b <- sprintf(paste0(
  "Scatterplot B compares the 7-dataset integrated dataset without Mikheev ",
  "(%d samples, %d DEGs) against the full 8-dataset integrated dataset ",
  "(%d samples, %d DEGs). ",
  "Of the %d DEGs in the 8-dataset integration, %d (%.0f%%) were also found ",
  "in the 7-dataset integration. ",
  "Among %d genes plotted (out of %d shared, at least one with FDR < 0.05), ",
  "%d were significant in both, %d in the 8-dataset integration only, ",
  "and %d in the 7-dataset integration only. ",
  "Log-fold-change agreement: Pearson r = %.3f, Spearman rho = %.3f, ",
  "same direction in %.1f%% of genes."),
  n_7m_filt, n_deg_7m, n_8ds_filt, n_deg_8ds,
  n_deg_8ds, n_7m_8ds_overlap, 100 * n_7m_8ds_overlap / max(n_deg_8ds, 1),
  nrow(df_b), n_all_b, n_both_b, n_8ds_b, n_7ds_b, r_b, rho_b, sd_b)

caption <- paste0(
  "Figure. Leave-one-out validation for the Mikheev dataset (GSE9984). ",
  "Only genes with FDR < 0.05 in at least one of the two compared analyses are shown. ",
  "The dotted line marks the identity line (y = x). ",
  "DEG thresholds: FDR < 0.05 and absolute log-fold-change greater than 1.\n\n",
  caption_a, "\n\n", caption_b)

docx_path <- file.path(output_dir, "scatter_mikheev_soncin_6_10w.docx")
doc <- read_docx()
doc <- body_add_img(doc, img_path, width = 7.5, height = 3.4)
doc <- body_add_par(doc, caption, style = "Normal")
print(doc, target = docx_path)
cat(sprintf("DOCX saved: %s\n", docx_path))
