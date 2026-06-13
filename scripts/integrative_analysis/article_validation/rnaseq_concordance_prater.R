options(bitmapType = "cairo")
library(ggplot2)
library(ggrepel)
library(openxlsx)

lin_ccc <- function(x, y) {
  ok <- complete.cases(x, y)
  if (sum(ok) < 3) return(NA_real_)
  x <- x[ok]; y <- y[ok]
  mx <- mean(x); my <- mean(y)
  sx <- var(x); sy <- var(y)
  sxy <- cov(x, y)
  2 * sxy / (sx + sy + (mx - my)^2)
}

output_dir <- "output/article_validation/rnaseq_concordance_prater"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

fig_dir <- "articles/imputation_article/figures"

PRATER_PATH <- "articles/yehor_conference_2026/references/Prater_2021_RNA-Seq_first-second_trimester_transition_supp_tables.xlsx"
FULL_BASE <- "output/yehor_sashko/phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko"
TRIM_BASE <- "output/yehor_sashko/phase2b_1_2_yehor_6ds_trimmed_to_prater"

# --- Load Prater 2021 RNA-seq DEGs ---

prater <- read.xlsx(loadWorkbook(PRATER_PATH),
                    sheet = "T1 DEGs_results_table_l2fc1")
prater$entrez <- as.character(prater$entrezgene_id)
prater <- prater[!is.na(prater$entrez) & prater$entrez != "NA", ]

prater_logfc <- setNames(as.numeric(prater$log2FoldChange), prater$entrez)
prater_padj  <- setNames(as.numeric(prater$padj), prater$entrez)
prater_sig   <- prater$entrez[!is.na(prater_padj) & prater_padj < 0.05 &
                               abs(prater_logfc[prater$entrez]) > 1]

cat(sprintf("Prater 2021: %d DEGs total, %d with ENTREZ ID, %d significant\n",
            nrow(read.xlsx(loadWorkbook(PRATER_PATH),
                           sheet = "T1 DEGs_results_table_l2fc1")),
            length(prater_logfc), length(prater_sig)))

# --- Helper: compute concordance stats for a pipeline run ---

compute_concordance <- function(de_path, sig_path, label) {
  full_de <- read.delim(de_path, stringsAsFactors = FALSE)
  sig_de  <- read.delim(sig_path, stringsAsFactors = FALSE)

  ours_logfc   <- setNames(full_de$logFC, full_de$gene)
  ours_fdr     <- setNames(full_de$adj.P.Val, full_de$gene)
  ours_sig_set <- sig_de$gene[sig_de$adj.P.Val < 0.05 & abs(sig_de$logFC) > 1]

  shared <- intersect(names(ours_logfc), names(prater_logfc))

  df <- data.frame(
    gene         = shared,
    ours_logfc   = ours_logfc[shared],
    prater_logfc = prater_logfc[shared],
    stringsAsFactors = FALSE
  )
  df$ours_sig   <- df$gene %in% ours_sig_set
  df$prater_sig <- df$gene %in% prater_sig
  df$same_dir   <- sign(df$ours_logfc) == sign(df$prater_logfc)

  r_val   <- cor(df$ours_logfc, df$prater_logfc, use = "complete.obs")
  ccc_val <- lin_ccc(df$ours_logfc, df$prater_logfc)

  both_sig <- sum(df$ours_sig & df$prater_sig)
  ours_only <- sum(df$ours_sig & !df$prater_sig)
  prater_only <- sum(!df$ours_sig & df$prater_sig)
  neither <- sum(!df$ours_sig & !df$prater_sig)

  both_sig_concordant <- sum(df$ours_sig & df$prater_sig & df$same_dir)
  all_concordant <- sum(df$same_dir)

  cat(sprintf("\n=== %s ===\n", label))
  cat(sprintf("Shared genes: %d\n", nrow(df)))
  cat(sprintf("Pearson r (logFC): %.3f\n", r_val))
  cat(sprintf("Lin's CCC (logFC): %.3f\n", ccc_val))
  cat(sprintf("Both significant: %d\n", both_sig))
  cat(sprintf("  directionally concordant: %d/%d (%.1f%%)\n",
              both_sig_concordant, both_sig,
              if (both_sig > 0) 100 * both_sig_concordant / both_sig else 0))
  cat(sprintf("Ours only: %d, Prater only: %d, Neither: %d\n",
              ours_only, prater_only, neither))
  cat(sprintf("All shared directionally concordant: %d/%d (%.1f%%)\n",
              all_concordant, nrow(df), 100 * all_concordant / nrow(df)))

  list(df = df, r = r_val, ccc = ccc_val, n_shared = nrow(df),
       both_sig = both_sig, both_sig_concordant = both_sig_concordant,
       ours_only = ours_only, prater_only = prater_only, neither = neither,
       all_concordant = all_concordant, label = label,
       n_testable = nrow(full_de), n_degs = nrow(sig_de))
}

# --- Compute for full and trimmed pipelines ---

full <- compute_concordance(
  file.path(FULL_BASE, "difexp_softimpute_combat_ref.tsv"),
  file.path(FULL_BASE, "difexp_significant_softimpute_combat_ref.tsv"),
  "Full pipeline (117 samples, 4-12 wk vs 14-19 wk)")

trim <- compute_concordance(
  file.path(TRIM_BASE, "difexp_softimpute_combat_ref.tsv"),
  file.path(TRIM_BASE, "difexp_significant_softimpute_combat_ref.tsv"),
  "Trimmed pipeline (36 samples, 7-8 wk vs 14 wk)")

# --- Scatterplot (full pipeline, article figure) ---

df <- full$df
df$category <- "Not significant in either"
df$category[df$ours_sig & df$prater_sig]  <- "Both significant"
df$category[df$ours_sig & !df$prater_sig] <- "This study only"
df$category[!df$ours_sig & df$prater_sig] <- "Prater 2021 only"

df$category <- factor(df$category,
  levels = c("Both significant", "This study only",
             "Prater 2021 only", "Not significant in either"))

cat_colors <- c("Both significant" = "#D32F2F",
                "This study only" = "#1976D2",
                "Prater 2021 only" = "#F57C00",
                "Not significant in either" = "#BDBDBD")

p <- ggplot(df, aes(x = ours_logfc, y = prater_logfc, color = category)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey40", linewidth = 0.4) +
  geom_hline(yintercept = c(-1, 1), linetype = "dotted",
             color = "steelblue", alpha = 0.4, linewidth = 0.3) +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted",
             color = "steelblue", alpha = 0.4, linewidth = 0.3) +
  geom_point(alpha = 0.5, size = 0.8) +
  scale_color_manual(values = cat_colors,
                     name = expression("Significance (FDR < 0.05, |logFC| > 1)")) +
  annotate("text", x = -3.5, y = 4.5,
           label = sprintf("r = %.3f\nCCC = %.3f\nn = %s",
                           full$r, full$ccc,
                           format(full$n_shared, big.mark = ",")),
           hjust = 0, vjust = 1, size = 3.5) +
  labs(x = "logFC (this study, microarray, 117 samples)",
       y = expression("log"[2]*"FC (Prater 2021, RNA-seq, 14 samples)")) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.key.size = unit(3, "mm")) +
  guides(color = guide_legend(nrow = 2, override.aes = list(size = 2, alpha = 0.8))) +
  scale_x_continuous(breaks = -4:4) +
  scale_y_continuous(breaks = -5:5) +
  coord_cartesian(xlim = c(-4, 4), ylim = c(-5, 5))

fig_path <- file.path(fig_dir, "fig_rnaseq_concordance_prater.png")
png(fig_path, width = 7, height = 6.5, units = "in", res = 300, type = "cairo")
print(p)
invisible(dev.off())
cat("\nSaved figure:", fig_path, "\n")

# --- Scatterplot (trimmed pipeline, article figure) ---

df_trim <- trim$df
df_trim$category <- "Not significant in either"
df_trim$category[df_trim$ours_sig & df_trim$prater_sig]  <- "Both significant"
df_trim$category[df_trim$ours_sig & !df_trim$prater_sig] <- "This study only"
df_trim$category[!df_trim$ours_sig & df_trim$prater_sig] <- "Prater 2021 only"

df_trim$category <- factor(df_trim$category,
  levels = c("Both significant", "This study only",
             "Prater 2021 only", "Not significant in either"))

p_trim <- ggplot(df_trim, aes(x = ours_logfc, y = prater_logfc,
                              color = category)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey60", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey40", linewidth = 0.4) +
  geom_hline(yintercept = c(-1, 1), linetype = "dotted",
             color = "steelblue", alpha = 0.4, linewidth = 0.3) +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted",
             color = "steelblue", alpha = 0.4, linewidth = 0.3) +
  geom_point(alpha = 0.5, size = 0.8) +
  scale_color_manual(values = cat_colors,
                     name = expression(
                       "Significance (FDR < 0.05, |logFC| > 1)")) +
  annotate("text", x = -3.5, y = 4.5,
           label = sprintf("r = %.3f\nCCC = %.3f\nn = %s",
                           trim$r, trim$ccc,
                           format(trim$n_shared, big.mark = ",")),
           hjust = 0, vjust = 1, size = 3.5) +
  labs(x = "logFC (this study, microarray, 36 samples, 7–8 wk)",
       y = expression(
         "log"[2]*"FC (Prater 2021, RNA-seq, 14 samples)")) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.key.size = unit(3, "mm")) +
  guides(color = guide_legend(nrow = 2,
         override.aes = list(size = 2, alpha = 0.8))) +
  scale_x_continuous(breaks = -4:4) +
  scale_y_continuous(breaks = -5:5) +
  coord_cartesian(xlim = c(-4, 4), ylim = c(-5, 5))

fig_trim_path <- file.path(fig_dir,
                           "fig_rnaseq_concordance_prater_trimmed.png")
png(fig_trim_path, width = 7, height = 6.5, units = "in",
    res = 300, type = "cairo")
print(p_trim)
invisible(dev.off())
cat("Saved trimmed figure:", fig_trim_path, "\n")

# --- Write summary for LaTeX ---

summary_text <- sprintf("
=== RNA-seq cross-validation summary (for LaTeX text) ===

Prater 2021: RNA-seq, 7-8 wk vs 13-14 wk, n=14, Biology Open 10(6):bio058222
  Total DEGs in Table T1: 3,268 (pre-filtered: padj<0.05, |log2FC|>1)
  With ENTREZ ID: %d

--- Full pipeline (117 samples, 4-12 wk 1T vs 14-19 wk 2T) ---
  Testable genes: %s
  DEGs (softImpute + ComBat-ref, FDR<0.05, |logFC|>1): %d
  Shared genes with Prater: %s
  Pearson r (logFC): %.3f
  Both significant (FDR<0.05, |logFC|>1 in both): %d
    directionally concordant: %d/%d (%.1f%%)
  This study only: %d
  Prater only: %d
  Neither significant: %d
  All shared genes directionally concordant: %d/%d (%.1f%%)

--- Trimmed pipeline (36 samples, 7-8 wk 1T vs 14 wk 2T) ---
  Testable genes: %s
  DEGs (softImpute + ComBat-ref, FDR<0.05, |logFC|>1): %d
  Shared genes with Prater: %s
  Pearson r (logFC): %.3f
  Both significant: %d
    directionally concordant: %d/%d (%.1f%%)
  This study only: %d
  Prater only: %d
  Neither significant: %d
  All shared genes directionally concordant: %d/%d (%.1f%%)
",
length(prater_logfc),
format(full$n_testable, big.mark = ","), full$n_degs,
format(full$n_shared, big.mark = ","), full$r,
full$both_sig, full$both_sig_concordant, full$both_sig,
if (full$both_sig > 0) 100 * full$both_sig_concordant / full$both_sig else 0,
full$ours_only, full$prater_only, full$neither,
full$all_concordant, full$n_shared, 100 * full$all_concordant / full$n_shared,
format(trim$n_testable, big.mark = ","), trim$n_degs,
format(trim$n_shared, big.mark = ","), trim$r,
trim$both_sig, trim$both_sig_concordant, trim$both_sig,
if (trim$both_sig > 0) 100 * trim$both_sig_concordant / trim$both_sig else 0,
trim$ours_only, trim$prater_only, trim$neither,
trim$all_concordant, trim$n_shared, 100 * trim$all_concordant / trim$n_shared
)

cat(summary_text)
writeLines(summary_text, file.path(output_dir, "concordance_summary.txt"))

write.csv(full$df, file.path(output_dir, "full_pipeline_vs_prater.csv"),
          row.names = FALSE)
write.csv(trim$df, file.path(output_dir, "trimmed_pipeline_vs_prater.csv"),
          row.names = FALSE)

cat("\nAll outputs saved to:", output_dir, "\n")
