options(bitmapType = "cairo")
library(ggplot2)
library(ggrepel)

output_dir <- "output/article_validation/prior_study_concordance"
base_dir <- "output/yehor_sashko/phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko"

# --- Load data ---

sig_softimpute <- read.delim(file.path(base_dir, "difexp_significant_softimpute_combat_ref.tsv"),
                             stringsAsFactors = FALSE)
sig_intersection <- read.delim(file.path(base_dir, "difexp_significant_none_combat_ref.tsv"),
                               stringsAsFactors = FALSE)
full_7ds <- read.delim(file.path(base_dir, "difexp_softimpute_combat_ref.tsv"),
                       stringsAsFactors = FALSE)
lykhenko_all <- read.csv("output/dissertation/integrative_pipeline_disser/1_2/difexp/difexp_all.csv",
                         stringsAsFactors = FALSE)
lykhenko_sig <- read.csv("output/dissertation/integrative_pipeline_disser/1_2/difexp/difexp_filtered.csv",
                         stringsAsFactors = FALSE)

lykhenko_all$ENTREZID <- as.character(lykhenko_all$ENTREZID)
lykhenko_sig$ENTREZID <- as.character(lykhenko_sig$ENTREZID)
full_7ds$gene <- as.character(full_7ds$gene)

# --- All 7ds DEGs (union of gained + intersection) ---

all_7ds_deg_ids <- unique(c(as.character(sig_softimpute$gene),
                            as.character(sig_intersection$gene)))
lykhenko_sig_ids <- lykhenko_sig$ENTREZID

# --- Classify every DEG from either analysis ---

all_deg_ids <- unique(c(all_7ds_deg_ids, lykhenko_sig_ids))

in_7ds     <- all_deg_ids %in% all_7ds_deg_ids
in_2021    <- all_deg_ids %in% lykhenko_sig_ids
in_2021_full <- all_deg_ids %in% lykhenko_all$ENTREZID
in_7ds_full  <- all_deg_ids %in% full_7ds$gene

# Three groups:
# 1. Both: DEG in 7ds AND DEG in 2021 → arrows
# 2. 7ds only: DEG in 7ds but NOT DEG in 2021 → single dot at 7ds position
# 3. 2021 only: DEG in 2021 but NOT DEG in 7ds → single dot at 2021 position

both_ids     <- all_deg_ids[in_7ds & in_2021]
only_7ds_ids <- all_deg_ids[in_7ds & !in_2021]
only_2021_ids <- all_deg_ids[!in_7ds & in_2021]

cat("DEGs in both analyses (arrows):", length(both_ids), "\n")
cat("DEGs only in 7ds:", length(only_7ds_ids), "\n")
cat("DEGs only in 2021:", length(only_2021_ids), "\n\n")

# --- Arrow data: genes DEG in both ---

# For 7ds position: use softImpute full table if available, else intersection full table
get_7ds_pos <- function(ids) {
  idx <- match(ids, full_7ds$gene)
  found <- !is.na(idx)
  logFC <- rep(NA_real_, length(ids))
  fdr   <- rep(NA_real_, length(ids))
  logFC[found] <- full_7ds$logFC[idx[found]]
  fdr[found]   <- full_7ds$adj.P.Val[idx[found]]
  list(logFC = logFC, fdr = fdr, found = found)
}

# Both: need 2021 and 7ds positions
both_in_2021 <- both_ids[both_ids %in% lykhenko_all$ENTREZID]
both_in_7ds  <- both_in_2021[both_in_2021 %in% full_7ds$gene]

idx_2021 <- match(both_in_7ds, lykhenko_all$ENTREZID)
idx_7ds  <- match(both_in_7ds, full_7ds$gene)

arrows_df <- data.frame(
  entrezid   = both_in_7ds,
  symbol     = lykhenko_all$SYMBOL[idx_2021],
  logFC_2021 = lykhenko_all$logFC[idx_2021],
  fdr_2021   = lykhenko_all$adj.P.Val[idx_2021],
  logFC_7ds  = full_7ds$logFC[idx_7ds],
  fdr_7ds    = full_7ds$adj.P.Val[idx_7ds],
  stringsAsFactors = FALSE
)

# --- 7ds-only dots ---

only_7ds_in_full <- only_7ds_ids[only_7ds_ids %in% full_7ds$gene]
idx_7ds_only <- match(only_7ds_in_full, full_7ds$gene)

dots_7ds_only <- data.frame(
  entrezid = only_7ds_in_full,
  logFC    = full_7ds$logFC[idx_7ds_only],
  fdr      = full_7ds$adj.P.Val[idx_7ds_only],
  stringsAsFactors = FALSE
)

# Also check which of these have a 2021 position (for possible labeling)
idx_in_2021 <- match(only_7ds_in_full, lykhenko_all$ENTREZID)
dots_7ds_only$symbol <- ifelse(!is.na(idx_in_2021),
                               lykhenko_all$SYMBOL[idx_in_2021], only_7ds_in_full)

# --- 2021-only dots ---

idx_2021_only <- match(only_2021_ids, lykhenko_all$ENTREZID)

dots_2021_only <- data.frame(
  entrezid = only_2021_ids,
  symbol   = lykhenko_all$SYMBOL[idx_2021_only],
  logFC    = lykhenko_all$logFC[idx_2021_only],
  fdr      = lykhenko_all$adj.P.Val[idx_2021_only],
  stringsAsFactors = FALSE
)

# --- Background: all 2021 genes ---

bg_df <- data.frame(
  logFC = lykhenko_all$logFC,
  nlog10_fdr = -log10(pmax(lykhenko_all$adj.P.Val, 1e-50)),
  stringsAsFactors = FALSE
)

# --- Cap FDR ---

fdr_cap <- 14
bg_df$nlog10_fdr <- pmin(bg_df$nlog10_fdr, fdr_cap)

arrows_df$nlog10_2021 <- pmin(-log10(pmax(arrows_df$fdr_2021, 1e-50)), fdr_cap)
arrows_df$nlog10_7ds  <- pmin(-log10(pmax(arrows_df$fdr_7ds, 1e-50)), fdr_cap)

dots_7ds_only$nlog10_fdr <- pmin(-log10(pmax(dots_7ds_only$fdr, 1e-50)), fdr_cap)
dots_2021_only$nlog10_fdr <- pmin(-log10(pmax(dots_2021_only$fdr, 1e-50)), fdr_cap)

# --- Labels ---

arrows_df$delta <- sqrt((arrows_df$logFC_7ds - arrows_df$logFC_2021)^2 +
                        (arrows_df$nlog10_7ds - arrows_df$nlog10_2021)^2)
top_arrows <- arrows_df[order(-arrows_df$delta), ][1:8, ]

top_7ds_only <- dots_7ds_only[order(-dots_7ds_only$nlog10_fdr), ][1:5, ]
top_2021_only <- dots_2021_only[order(-dots_2021_only$nlog10_fdr), ][1:5, ]

# --- Legend labels with counts ---

arrow_label   <- sprintf("DEG in both (n=%d)", nrow(arrows_df))
only7ds_label <- sprintf("DEG only in 7ds (n=%d)", nrow(dots_7ds_only))
only2021_label <- sprintf("DEG only in 2021 (n=%d)", nrow(dots_2021_only))

colours <- c("#4393C3", "#D6604D", "#5AAE61")
names(colours) <- c(arrow_label, only7ds_label, only2021_label)

# --- Plot ---

p <- ggplot() +
  # background
  geom_point(data = bg_df, aes(x = logFC, y = nlog10_fdr),
             colour = "grey88", size = 0.3, alpha = 0.2) +
  # threshold lines
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey50") +
  # arrows: DEG in both analyses (one colour)
  geom_segment(data = arrows_df,
               aes(x = logFC_2021, y = nlog10_2021,
                   xend = logFC_7ds, yend = nlog10_7ds,
                   colour = arrow_label),
               arrow = arrow(length = unit(1.2, "mm"), type = "closed"),
               alpha = 0.25, linewidth = 0.35) +
  # 2021 end of arrows (open)
  geom_point(data = arrows_df,
             aes(x = logFC_2021, y = nlog10_2021, colour = arrow_label),
             shape = 1, size = 1.2, alpha = 0.4) +
  # 7ds end of arrows (filled)
  geom_point(data = arrows_df,
             aes(x = logFC_7ds, y = nlog10_7ds, colour = arrow_label),
             shape = 16, size = 1.6, alpha = 0.6) +
  # 7ds-only DEGs (no arrow, just dot)
  geom_point(data = dots_7ds_only,
             aes(x = logFC, y = nlog10_fdr, colour = only7ds_label),
             shape = 17, size = 2.2, alpha = 0.7) +
  # 2021-only DEGs (no arrow, just dot)
  geom_point(data = dots_2021_only,
             aes(x = logFC, y = nlog10_fdr, colour = only2021_label),
             shape = 15, size = 2.2, alpha = 0.7) +
  # labels: top arrows
  geom_text_repel(data = top_arrows,
                  aes(x = logFC_7ds, y = nlog10_7ds, label = symbol),
                  size = 2.8, fontface = "italic", colour = colours[1],
                  max.overlaps = 20, min.segment.length = 0) +
  # labels: top 7ds-only
  geom_text_repel(data = top_7ds_only,
                  aes(x = logFC, y = nlog10_fdr, label = symbol),
                  size = 2.8, fontface = "italic", colour = colours[2],
                  max.overlaps = 20, min.segment.length = 0) +
  # labels: top 2021-only
  geom_text_repel(data = top_2021_only,
                  aes(x = logFC, y = nlog10_fdr, label = symbol),
                  size = 2.8, fontface = "italic", colour = colours[3],
                  max.overlaps = 20, min.segment.length = 0) +
  scale_colour_manual(values = colours, name = NULL) +
  labs(x = expression(log[2]~"fold change"),
       y = expression(-log[10]~"FDR"),
       title = "DEG transitions between 2021 (4 Affy, 22 samples) and 7ds integration (123 samples)",
       subtitle = expression("Arrows: DEG in both analyses (open circle"
                             %->% "filled circle). Triangles/squares: DEG in one analysis only.")) +
  coord_cartesian(xlim = c(-5.5, 4.5), ylim = c(0, fdr_cap + 0.5)) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 11)) +
  guides(colour = guide_legend(ncol = 1,
    override.aes = list(size = 3, alpha = 1, linewidth = 0.8,
                        shape = c(16, 17, 15))))

ggsave(file.path(output_dir, "fig_volcano_arrows_simple.pdf"), p, width = 11, height = 9)
png(file.path(output_dir, "fig_volcano_arrows_simple.png"),
    width = 3300, height = 2700, res = 300, type = "cairo")
print(p)
invisible(dev.off())

cat("\nSaved to:", file.path(output_dir, "fig_volcano_arrows_simple.{pdf,png}"), "\n")
