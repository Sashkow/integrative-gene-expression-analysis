#!/usr/bin/env Rscript
# Side-by-side scatters for Soncin with 1T restricted to weeks 6-10
# A: GSE100051 single vs 8ds_soncin_6_10w
# B: 8ds_soncin_6_10w vs 8ds original (full Soncin)

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

deg_soncin_unfilt <- read.csv(
  "output/volodymyr/volodymyr_1_2_trim/first_GSE100051_single_difexp_unfiltered.csv",
  stringsAsFactors = FALSE)
deg_6_10w_full <- read.delim(
  "output/yehor_sashko/phase2b_1_2_yehor_8ds_soncin_6_10w/difexp_softimpute_combat_ref.tsv",
  stringsAsFactors = FALSE)
deg_8ds_full <- read.delim(
  "output/yehor_sashko/phase2b_1_2_yehor_8ds_enriched_sashko/difexp_softimpute_combat_ref.tsv",
  stringsAsFactors = FALSE)

# symbol lookup
id2sym <- setNames(deg_soncin_unfilt$SYMBOL, deg_soncin_unfilt$ENTREZID)

# sample counts
pheno <- read.delim(
  "data/mapped/yehor/8_yehor_preprocessed_datasets_enriched_sashko/phenodata_placenta_1_2_enriched_sashko.tsv",
  stringsAsFactors = FALSE)
pheno <- pheno[pheno$condition == "healthy" &
  pheno$trimester %in% c("First trimester", "Second trimester"), ]
n_total <- nrow(pheno)
son_pheno <- pheno[pheno$dataset_id == "GSE100051", ]
n_son_orig <- nrow(son_pheno)
n_son_6_10w <- sum(
  (son_pheno$trimester == "First trimester" &
   son_pheno$fetus_week >= 6 & son_pheno$fetus_week <= 10) |
  son_pheno$trimester == "Second trimester")
n_filtered <- n_total - n_son_orig + n_son_6_10w
n_dropped <- n_son_orig - n_son_6_10w

# --- Plot A: GSE100051 single vs 8ds_soncin_6_10w ---
merged_a <- merge(
  deg_soncin_unfilt[, c("ENTREZID", "SYMBOL", "logFC", "P.Value", "adj.P.Val")],
  deg_6_10w_full[, c("gene", "logFC", "P.Value", "adj.P.Val")],
  by.x = "ENTREZID", by.y = "gene", suffixes = c("_single", "_6_10w"))

n_all_a <- nrow(merged_a)
df_a <- merged_a[merged_a$adj.P.Val_single < 0.05 |
                 merged_a$adj.P.Val_6_10w < 0.05, ]

df_a$sig_cat <- ifelse(
  df_a$adj.P.Val_single < 0.05 & df_a$adj.P.Val_6_10w < 0.05,
  "both FDR < 0.05",
  ifelse(df_a$adj.P.Val_6_10w < 0.05, "8ds filtered only",
         "GSE100051 only"))

r_a   <- cor(df_a$logFC_single, df_a$logFC_6_10w, method = "pearson")
rho_a <- cor(df_a$logFC_single, df_a$logFC_6_10w, method = "spearman")
sd_a  <- mean(sign(df_a$logFC_single) == sign(df_a$logFC_6_10w) &
              df_a$logFC_single != 0 & df_a$logFC_6_10w != 0) * 100
n_both_a <- sum(df_a$sig_cat == "both FDR < 0.05")
n_filt_a <- sum(df_a$sig_cat == "8ds filtered only")
n_sing_a <- sum(df_a$sig_cat == "GSE100051 only")

df_a$dist <- sqrt(df_a$logFC_single^2 + df_a$logFC_6_10w^2)
df_a$label <- ifelse(rank(-df_a$dist) <= 15, df_a$SYMBOL, NA)

# --- Plot B: 8ds_soncin_6_10w vs 8ds original ---
merged_b <- merge(
  deg_6_10w_full[, c("gene", "logFC", "P.Value", "adj.P.Val")],
  deg_8ds_full[, c("gene", "logFC", "P.Value", "adj.P.Val")],
  by = "gene", suffixes = c("_6_10w", "_8ds"))
merged_b$SYMBOL <- id2sym[as.character(merged_b$gene)]

n_all_b <- nrow(merged_b)
df_b <- merged_b[merged_b$adj.P.Val_6_10w < 0.05 |
                 merged_b$adj.P.Val_8ds < 0.05, ]

df_b$sig_cat <- ifelse(
  df_b$adj.P.Val_6_10w < 0.05 & df_b$adj.P.Val_8ds < 0.05,
  "both FDR < 0.05",
  ifelse(df_b$adj.P.Val_8ds < 0.05, "8ds original only",
         "8ds filtered only"))

r_b   <- cor(df_b$logFC_6_10w, df_b$logFC_8ds, method = "pearson")
rho_b <- cor(df_b$logFC_6_10w, df_b$logFC_8ds, method = "spearman")
sd_b  <- mean(sign(df_b$logFC_6_10w) == sign(df_b$logFC_8ds) &
              df_b$logFC_6_10w != 0 & df_b$logFC_8ds != 0) * 100
n_both_b <- sum(df_b$sig_cat == "both FDR < 0.05")
n_orig_b <- sum(df_b$sig_cat == "8ds original only")
n_filt_b <- sum(df_b$sig_cat == "8ds filtered only")

df_b$dist <- sqrt(df_b$logFC_6_10w^2 + df_b$logFC_8ds^2)
df_b$label <- ifelse(rank(-df_b$dist) <= 15, df_b$SYMBOL, NA)

# ============================================================
# Plots
# ============================================================

sig_colors_a <- c("both FDR < 0.05" = "#E69F00",
                  "8ds filtered only" = "#56B4E9",
                  "GSE100051 only" = "#009E73")

pa <- ggplot(df_a, aes(x = logFC_single, y = logFC_6_10w, color = sig_cat)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey40") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20,
                  show.legend = FALSE) +
  scale_color_manual(values = sig_colors_a, name = NULL) +
  labs(title = "A. GSE100051 vs 8ds (Soncin 6-10w)",
       x = sprintf("GSE100051 logFC (%d samples)", n_son_orig),
       y = sprintf("8ds logFC, Soncin 6-10w (%d samples)", n_filtered)) +
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
                  "8ds original only" = "#009E73",
                  "8ds filtered only" = "#56B4E9")

pb <- ggplot(df_b, aes(x = logFC_6_10w, y = logFC_8ds, color = sig_cat)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey40") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20,
                  show.legend = FALSE) +
  scale_color_manual(values = sig_colors_b, name = NULL) +
  labs(title = "B. 8ds (Soncin 6-10w) vs 8ds original",
       x = sprintf("8ds logFC, Soncin 6-10w (%d samples)", n_filtered),
       y = sprintf("8ds logFC, all Soncin (%d samples)", n_total)) +
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

img_path <- file.path(output_dir, "scatter_soncin_6_10w_side_by_side.png")
png(img_path, width = 2400, height = 1100, res = 150)
print(combined)
dev.off()

# ============================================================
# DOCX
# ============================================================

deg_6_10w_sig <- read.delim(
  "output/yehor_sashko/phase2b_1_2_yehor_8ds_soncin_6_10w/difexp_significant_softimpute_combat_ref.tsv",
  stringsAsFactors = FALSE)
deg_8ds_sig <- read.delim(
  "output/yehor_sashko/phase2b_1_2_yehor_8ds_enriched_sashko/difexp_significant_softimpute_combat_ref.tsv",
  stringsAsFactors = FALSE)
n_deg_6_10w <- nrow(deg_6_10w_sig)
n_deg_8ds <- nrow(deg_8ds_sig)
n_deg_overlap <- length(intersect(
  as.character(deg_6_10w_sig$gene), as.character(deg_8ds_sig$gene)))

caption_a <- sprintf(paste0(
  "A. GSE100051 single-dataset (limma, %d samples) vs 8ds with Soncin ",
  "restricted to gestational weeks 6-10 (%d samples, %d Soncin 1T samples dropped). ",
  "Genes plotted: %d / %d (at least one FDR < 0.05). ",
  "Both sig.: %d, 8ds filtered only: %d, GSE100051 only: %d. ",
  "Pearson r = %.3f, Spearman rho = %.3f, same direction = %.1f%%."),
  n_son_orig, n_filtered, n_dropped,
  nrow(df_a), n_all_a, n_both_a, n_filt_a, n_sing_a, r_a, rho_a, sd_a)

caption_b <- sprintf(paste0(
  "B. 8ds with Soncin weeks 6-10 (%d samples, %d DEGs) vs 8ds original ",
  "(%d samples, %d DEGs). DEG overlap: %d (%.0f%% of original retained). ",
  "Genes plotted: %d / %d (at least one FDR < 0.05). ",
  "Both sig.: %d, 8ds original only: %d, 8ds filtered only: %d. ",
  "Pearson r = %.3f, Spearman rho = %.3f, same direction = %.1f%%."),
  n_filtered, n_deg_6_10w, n_total, n_deg_8ds,
  n_deg_overlap, 100 * n_deg_overlap / n_deg_8ds,
  nrow(df_b), n_all_b, n_both_b, n_orig_b, n_filt_b, r_b, rho_b, sd_b)

caption <- paste0(
  "Figure. Effect of restricting Soncin (GSE100051) first-trimester samples ",
  "to gestational weeks 6-10 (removing weeks 4, 5, 11, 12). ",
  "Only genes with FDR < 0.05 in at least one comparison are shown. ",
  "Dotted line = identity (y = x). DEG thresholds: FDR < 0.05 and |logFC| > 1. ",
  "Pipeline: softimpute + combat_ref.\n\n",
  caption_a, "\n\n", caption_b)

docx_path <- file.path(output_dir, "scatter_soncin_6_10w.docx")
doc <- read_docx()
doc <- body_add_img(doc, img_path, width = 7.5, height = 3.4)
doc <- body_add_par(doc, caption, style = "Normal")
print(doc, target = docx_path)
cat(sprintf("DOCX saved: %s\n", docx_path))
