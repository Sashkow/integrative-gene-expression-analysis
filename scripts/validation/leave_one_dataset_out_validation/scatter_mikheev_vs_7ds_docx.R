#!/usr/bin/env Rscript
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
  "output/yehor_sashko/phase2b_1_2_yehor_7ds_no_mikheev/difexp_softimpute_combat_ref.tsv",
  stringsAsFactors = FALSE)
deg_8ds_full <- read.delim(
  "output/yehor_sashko/phase2b_1_2_yehor_8ds_enriched_sashko/difexp_softimpute_combat_ref.tsv",
  stringsAsFactors = FALSE)

# --- Plot 1 data: GSE9984 single vs 7ds no Mikheev ---
merged <- merge(
  deg_mikheev_unfilt[, c("ENTREZID", "SYMBOL", "logFC", "P.Value", "adj.P.Val")],
  deg_7m_full[, c("gene", "logFC", "P.Value", "adj.P.Val")],
  by.x = "ENTREZID", by.y = "gene", suffixes = c("_single", "_7ds"))

n_all <- nrow(merged)
df <- merged[merged$adj.P.Val_single < 0.05 | merged$adj.P.Val_7ds < 0.05, ]

df$sig_cat <- ifelse(
  df$adj.P.Val_single < 0.05 & df$adj.P.Val_7ds < 0.05,
  "both FDR < 0.05",
  ifelse(df$adj.P.Val_7ds < 0.05, "7ds only", "GSE9984 only"))

r1   <- cor(df$logFC_single, df$logFC_7ds, method = "pearson")
rho1 <- cor(df$logFC_single, df$logFC_7ds, method = "spearman")
sd1  <- mean(sign(df$logFC_single) == sign(df$logFC_7ds) &
             df$logFC_single != 0 & df$logFC_7ds != 0) * 100
n_both1  <- sum(df$sig_cat == "both FDR < 0.05")
n_7ds1   <- sum(df$sig_cat == "7ds only")
n_single <- sum(df$sig_cat == "GSE9984 only")

df$dist <- sqrt(df$logFC_single^2 + df$logFC_7ds^2)
df$label <- ifelse(rank(-df$dist) <= 15, df$SYMBOL, NA)

# --- Plot 2 data: 7ds no Mikheev vs 8ds ---
merged2 <- merge(
  deg_7m_full[, c("gene", "logFC", "P.Value", "adj.P.Val")],
  deg_8ds_full[, c("gene", "logFC", "P.Value", "adj.P.Val")],
  by = "gene", suffixes = c("_7ds", "_8ds"))
merged2$SYMBOL <- deg_mikheev_unfilt$SYMBOL[
  match(merged2$gene, deg_mikheev_unfilt$ENTREZID)]

n_all2 <- nrow(merged2)
df2 <- merged2[merged2$adj.P.Val_7ds < 0.05 | merged2$adj.P.Val_8ds < 0.05, ]

df2$sig_cat <- ifelse(
  df2$adj.P.Val_7ds < 0.05 & df2$adj.P.Val_8ds < 0.05,
  "both FDR < 0.05",
  ifelse(df2$adj.P.Val_8ds < 0.05, "8ds only", "7ds only"))

r2   <- cor(df2$logFC_7ds, df2$logFC_8ds, method = "pearson")
rho2 <- cor(df2$logFC_7ds, df2$logFC_8ds, method = "spearman")
sd2  <- mean(sign(df2$logFC_7ds) == sign(df2$logFC_8ds) &
             df2$logFC_7ds != 0 & df2$logFC_8ds != 0) * 100
n_both2 <- sum(df2$sig_cat == "both FDR < 0.05")
n_8ds   <- sum(df2$sig_cat == "8ds only")
n_7ds2  <- sum(df2$sig_cat == "7ds only")

df2$dist <- sqrt(df2$logFC_7ds^2 + df2$logFC_8ds^2)
df2$label <- ifelse(rank(-df2$dist) <= 15, df2$SYMBOL, NA)

# ============================================================
# Plots (minimal text on plot, stats go to caption)
# ============================================================

sig_colors <- c("both FDR < 0.05" = "#E69F00",
                "7ds only" = "#56B4E9",
                "GSE9984 only" = "#009E73")

p1 <- ggplot(df, aes(x = logFC_single, y = logFC_7ds, color = sig_cat)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey40") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20,
                  show.legend = FALSE) +
  scale_color_manual(values = sig_colors, name = NULL) +
  labs(title = "A. GSE9984 vs 7ds no Mikheev",
       x = "GSE9984 logFC (8 samples)",
       y = "7ds no Mikheev logFC (140 samples)") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 13),
        legend.position = "inside",
        legend.position.inside = c(0.22, 0.88),
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.35, "cm"),
        legend.background = element_rect(fill = alpha("white", 0.7),
                                         color = NA),
        plot.margin = margin(5, 10, 5, 5))

sig_colors2 <- c("both FDR < 0.05" = "#E69F00",
                 "8ds only" = "#009E73",
                 "7ds only" = "#56B4E9")

p2 <- ggplot(df2, aes(x = logFC_7ds, y = logFC_8ds, color = sig_cat)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey40") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20,
                  show.legend = FALSE) +
  scale_color_manual(values = sig_colors2, name = NULL) +
  labs(title = "B. 7ds no Mikheev vs 8ds",
       x = "7ds no Mikheev logFC (140 samples)",
       y = "8ds logFC (148 samples)") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 13),
        legend.position = "inside",
        legend.position.inside = c(0.22, 0.88),
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.35, "cm"),
        legend.background = element_rect(fill = alpha("white", 0.7),
                                         color = NA),
        plot.margin = margin(5, 5, 5, 10))

combined <- p1 + p2

img_path <- file.path(output_dir, "scatter_mikheev_7ds_8ds_side_by_side.png")
png(img_path, width = 2400, height = 1100, res = 150)
print(combined)
dev.off()

# ============================================================
# Soncin: GSE100051 single vs 7ds no Soncin + 7ds vs 8ds
# ============================================================

deg_soncin_unfilt <- read.csv(
  "output/volodymyr/volodymyr_1_2_trim/first_GSE100051_single_difexp_unfiltered.csv",
  stringsAsFactors = FALSE)
deg_7s_full <- read.delim(
  "output/yehor_sashko/phase2b_1_2_yehor_7ds_no_soncin/difexp_softimpute_combat_ref.tsv",
  stringsAsFactors = FALSE)

# symbol lookup combining both sources
id2sym <- setNames(deg_mikheev_unfilt$SYMBOL, deg_mikheev_unfilt$ENTREZID)
id2sym_s <- setNames(deg_soncin_unfilt$SYMBOL, deg_soncin_unfilt$ENTREZID)
id2sym[names(id2sym_s)] <- ifelse(
  is.na(id2sym[names(id2sym_s)]), id2sym_s, id2sym[names(id2sym_s)])

# --- Soncin plot A: GSE100051 single vs 7ds no Soncin ---
merged_s <- merge(
  deg_soncin_unfilt[, c("ENTREZID", "SYMBOL", "logFC", "P.Value", "adj.P.Val")],
  deg_7s_full[, c("gene", "logFC", "P.Value", "adj.P.Val")],
  by.x = "ENTREZID", by.y = "gene", suffixes = c("_single", "_7ds"))

n_all_s <- nrow(merged_s)
df_s <- merged_s[merged_s$adj.P.Val_single < 0.05 |
                 merged_s$adj.P.Val_7ds < 0.05, ]

df_s$sig_cat <- ifelse(
  df_s$adj.P.Val_single < 0.05 & df_s$adj.P.Val_7ds < 0.05,
  "both FDR < 0.05",
  ifelse(df_s$adj.P.Val_7ds < 0.05, "7ds only", "GSE100051 only"))

r_s1   <- cor(df_s$logFC_single, df_s$logFC_7ds, method = "pearson")
rho_s1 <- cor(df_s$logFC_single, df_s$logFC_7ds, method = "spearman")
sd_s1  <- mean(sign(df_s$logFC_single) == sign(df_s$logFC_7ds) &
               df_s$logFC_single != 0 & df_s$logFC_7ds != 0) * 100
n_both_s1  <- sum(df_s$sig_cat == "both FDR < 0.05")
n_7ds_s1   <- sum(df_s$sig_cat == "7ds only")
n_single_s <- sum(df_s$sig_cat == "GSE100051 only")

df_s$dist <- sqrt(df_s$logFC_single^2 + df_s$logFC_7ds^2)
df_s$label <- ifelse(rank(-df_s$dist) <= 15, df_s$SYMBOL, NA)

# --- Soncin plot B: 7ds no Soncin vs 8ds ---
merged_s2 <- merge(
  deg_7s_full[, c("gene", "logFC", "P.Value", "adj.P.Val")],
  deg_8ds_full[, c("gene", "logFC", "P.Value", "adj.P.Val")],
  by = "gene", suffixes = c("_7ds", "_8ds"))
merged_s2$SYMBOL <- id2sym[as.character(merged_s2$gene)]

n_all_s2 <- nrow(merged_s2)
df_s2 <- merged_s2[merged_s2$adj.P.Val_7ds < 0.05 |
                   merged_s2$adj.P.Val_8ds < 0.05, ]

df_s2$sig_cat <- ifelse(
  df_s2$adj.P.Val_7ds < 0.05 & df_s2$adj.P.Val_8ds < 0.05,
  "both FDR < 0.05",
  ifelse(df_s2$adj.P.Val_8ds < 0.05, "8ds only", "7ds only"))

r_s2   <- cor(df_s2$logFC_7ds, df_s2$logFC_8ds, method = "pearson")
rho_s2 <- cor(df_s2$logFC_7ds, df_s2$logFC_8ds, method = "spearman")
sd_s2  <- mean(sign(df_s2$logFC_7ds) == sign(df_s2$logFC_8ds) &
               df_s2$logFC_7ds != 0 & df_s2$logFC_8ds != 0) * 100
n_both_s2 <- sum(df_s2$sig_cat == "both FDR < 0.05")
n_8ds_s   <- sum(df_s2$sig_cat == "8ds only")
n_7ds_s2  <- sum(df_s2$sig_cat == "7ds only")

df_s2$dist <- sqrt(df_s2$logFC_7ds^2 + df_s2$logFC_8ds^2)
df_s2$label <- ifelse(rank(-df_s2$dist) <= 15, df_s2$SYMBOL, NA)

# --- Soncin plots ---
sig_colors_s <- c("both FDR < 0.05" = "#E69F00",
                  "7ds only" = "#56B4E9",
                  "GSE100051 only" = "#009E73")

ps1 <- ggplot(df_s, aes(x = logFC_single, y = logFC_7ds, color = sig_cat)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey40") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20,
                  show.legend = FALSE) +
  scale_color_manual(values = sig_colors_s, name = NULL) +
  labs(title = "A. GSE100051 vs 7ds no Soncin",
       x = "GSE100051 logFC (49 samples)",
       y = "7ds no Soncin logFC (99 samples)") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 13),
        legend.position = "inside",
        legend.position.inside = c(0.22, 0.88),
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.35, "cm"),
        legend.background = element_rect(fill = alpha("white", 0.7),
                                         color = NA),
        plot.margin = margin(5, 10, 5, 5))

sig_colors_s2 <- c("both FDR < 0.05" = "#E69F00",
                   "8ds only" = "#009E73",
                   "7ds only" = "#56B4E9")

ps2 <- ggplot(df_s2, aes(x = logFC_7ds, y = logFC_8ds, color = sig_cat)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey40") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20,
                  show.legend = FALSE) +
  scale_color_manual(values = sig_colors_s2, name = NULL) +
  labs(title = "B. 7ds no Soncin vs 8ds",
       x = "7ds no Soncin logFC (99 samples)",
       y = "8ds logFC (148 samples)") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 13),
        legend.position = "inside",
        legend.position.inside = c(0.22, 0.88),
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.35, "cm"),
        legend.background = element_rect(fill = alpha("white", 0.7),
                                         color = NA),
        plot.margin = margin(5, 5, 5, 10))

combined_s <- ps1 + ps2

img_path_s <- file.path(output_dir, "scatter_soncin_7ds_8ds_side_by_side.png")
png(img_path_s, width = 2400, height = 1100, res = 150)
print(combined_s)
dev.off()

# ============================================================
# DOCX
# ============================================================

caption_a <- sprintf(paste0(
  "A. GSE9984 single-dataset (limma, 8 samples) vs 7ds no Mikheev ",
  "(softimpute + combat_ref, 140 samples). ",
  "Genes plotted: %d / %d (at least one FDR < 0.05). ",
  "Both sig.: %d, 7ds only: %d, GSE9984 only: %d. ",
  "Pearson r = %.3f, Spearman rho = %.3f, same direction = %.1f%%."),
  nrow(df), n_all, n_both1, n_7ds1, n_single, r1, rho1, sd1)

caption_b <- sprintf(paste0(
  "B. 7ds no Mikheev (140 samples) vs 8ds full integration (148 samples), ",
  "both softimpute + combat_ref. ",
  "Genes plotted: %d / %d (at least one FDR < 0.05). ",
  "Both sig.: %d, 8ds only: %d, 7ds only: %d. ",
  "Pearson r = %.3f, Spearman rho = %.3f, same direction = %.1f%%."),
  nrow(df2), n_all2, n_both2, n_8ds, n_7ds2, r2, rho2, sd2)

caption_mik <- paste0(
  "Figure 1. Leave-one-out validation: effect of removing Mikheev (GSE9984) ",
  "on 1T vs 2T differential expression. ",
  "Only genes with FDR < 0.05 in at least one comparison are shown. ",
  "Dotted line = identity (y = x). DEG thresholds: FDR < 0.05 and |logFC| > 1.\n\n",
  caption_a, "\n\n", caption_b)

# --- Soncin captions ---
caption_sa <- sprintf(paste0(
  "A. GSE100051 single-dataset (limma, 49 samples) vs 7ds no Soncin ",
  "(softimpute + combat_ref, 99 samples). ",
  "Genes plotted: %d / %d (at least one FDR < 0.05). ",
  "Both sig.: %d, 7ds only: %d, GSE100051 only: %d. ",
  "Pearson r = %.3f, Spearman rho = %.3f, same direction = %.1f%%."),
  nrow(df_s), n_all_s, n_both_s1, n_7ds_s1, n_single_s, r_s1, rho_s1, sd_s1)

caption_sb <- sprintf(paste0(
  "B. 7ds no Soncin (99 samples) vs 8ds full integration (148 samples), ",
  "both softimpute + combat_ref. ",
  "Genes plotted: %d / %d (at least one FDR < 0.05). ",
  "Both sig.: %d, 8ds only: %d, 7ds only: %d. ",
  "Pearson r = %.3f, Spearman rho = %.3f, same direction = %.1f%%."),
  nrow(df_s2), n_all_s2, n_both_s2, n_8ds_s, n_7ds_s2, r_s2, rho_s2, sd_s2)

caption_son <- paste0(
  "Figure 2. Leave-one-out validation: effect of removing Soncin (GSE100051) ",
  "on 1T vs 2T differential expression. ",
  "Only genes with FDR < 0.05 in at least one comparison are shown. ",
  "Dotted line = identity (y = x). DEG thresholds: FDR < 0.05 and |logFC| > 1.\n\n",
  caption_sa, "\n\n", caption_sb)

docx_path <- file.path(output_dir, "scatter_leave_one_out.docx")
doc <- read_docx()
doc <- body_add_img(doc, img_path, width = 7.5, height = 3.4)
doc <- body_add_par(doc, caption_mik, style = "Normal")
doc <- body_add_break(doc)
doc <- body_add_img(doc, img_path_s, width = 7.5, height = 3.4)
doc <- body_add_par(doc, caption_son, style = "Normal")
print(doc, target = docx_path)
cat(sprintf("DOCX saved: %s\n", docx_path))
