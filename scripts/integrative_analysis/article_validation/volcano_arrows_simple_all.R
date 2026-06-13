options(bitmapType = "cairo")
library(ggplot2)
library(ggrepel)

output_dir <- "output/article_validation/prior_study_concordance"
base_dir <- "output/yehor_sashko/phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko"

# --- Load full limma tables ---

full_7ds <- read.delim(file.path(base_dir, "difexp_softimpute_combat_ref.tsv"),
                       stringsAsFactors = FALSE)
lykhenko_all <- read.csv("output/dissertation/integrative_pipeline_disser/1_2/difexp/difexp_all.csv",
                         stringsAsFactors = FALSE)

lykhenko_all$ENTREZID <- as.character(lykhenko_all$ENTREZID)
full_7ds$gene <- as.character(full_7ds$gene)

# --- DEG sets: FDR < 0.05 only (no logFC filter) ---

deg_7ds   <- full_7ds$gene[full_7ds$adj.P.Val < 0.05]
deg_2021  <- lykhenko_all$ENTREZID[lykhenko_all$adj.P.Val < 0.05]

all_degs <- unique(c(deg_7ds, deg_2021))

in_7ds  <- all_degs %in% deg_7ds
in_2021 <- all_degs %in% deg_2021

both_ids      <- all_degs[in_7ds & in_2021]
only_7ds_ids  <- all_degs[in_7ds & !in_2021]
only_2021_ids <- all_degs[!in_7ds & in_2021]

# Filter to genes present in both full tables (for position data)
both_ids      <- both_ids[both_ids %in% lykhenko_all$ENTREZID & both_ids %in% full_7ds$gene]
only_7ds_ids  <- only_7ds_ids[only_7ds_ids %in% full_7ds$gene]
only_2021_ids <- only_2021_ids[only_2021_ids %in% lykhenko_all$ENTREZID]

cat("DEGs in both (arrows):", length(both_ids), "\n")
cat("DEGs only in 7ds:", length(only_7ds_ids), "\n")
cat("DEGs only in 2021:", length(only_2021_ids), "\n\n")

# --- Arrow data ---

fdr_cap <- 14

idx_2021 <- match(both_ids, lykhenko_all$ENTREZID)
idx_7ds  <- match(both_ids, full_7ds$gene)

arrows_df <- data.frame(
  x0 = lykhenko_all$logFC[idx_2021],
  y0 = pmin(-log10(pmax(lykhenko_all$adj.P.Val[idx_2021], 1e-50)), fdr_cap),
  x1 = full_7ds$logFC[idx_7ds],
  y1 = pmin(-log10(pmax(full_7ds$adj.P.Val[idx_7ds], 1e-50)), fdr_cap),
  symbol = lykhenko_all$SYMBOL[idx_2021],
  stringsAsFactors = FALSE
)

# --- 7ds-only dots ---

idx_7ds_only <- match(only_7ds_ids, full_7ds$gene)
dots_7ds_only <- data.frame(
  logFC = full_7ds$logFC[idx_7ds_only],
  nlog10_fdr = pmin(-log10(pmax(full_7ds$adj.P.Val[idx_7ds_only], 1e-50)), fdr_cap),
  stringsAsFactors = FALSE
)

# Try to get symbols
idx_sym <- match(only_7ds_ids, lykhenko_all$ENTREZID)
dots_7ds_only$symbol <- ifelse(!is.na(idx_sym), lykhenko_all$SYMBOL[idx_sym], only_7ds_ids)

# --- 2021-only dots ---

idx_2021_only <- match(only_2021_ids, lykhenko_all$ENTREZID)
dots_2021_only <- data.frame(
  logFC = lykhenko_all$logFC[idx_2021_only],
  nlog10_fdr = pmin(-log10(pmax(lykhenko_all$adj.P.Val[idx_2021_only], 1e-50)), fdr_cap),
  symbol = lykhenko_all$SYMBOL[idx_2021_only],
  stringsAsFactors = FALSE
)

# --- Background: all 2021 genes ---

bg_df <- data.frame(
  logFC = lykhenko_all$logFC,
  nlog10_fdr = pmin(-log10(pmax(lykhenko_all$adj.P.Val, 1e-50)), fdr_cap),
  stringsAsFactors = FALSE
)

# --- Labels for top genes ---

arrows_df$delta <- sqrt((arrows_df$x1 - arrows_df$x0)^2 + (arrows_df$y1 - arrows_df$y0)^2)
top_arrows <- arrows_df[order(-arrows_df$delta), ][1:8, ]

top_7ds_only <- dots_7ds_only[order(-dots_7ds_only$nlog10_fdr), ][1:5, ]
top_2021_only <- dots_2021_only[order(-dots_2021_only$nlog10_fdr), ][1:5, ]

# --- Legend labels ---

arrow_label    <- sprintf("DEG in both (n=%d)", length(both_ids))
only7ds_label  <- sprintf("DEG only in 7ds (n=%d)", length(only_7ds_ids))
only2021_label <- sprintf("DEG only in 2021 (n=%d)", length(only_2021_ids))

colours <- c("#4393C3", "#D6604D", "#5AAE61")
names(colours) <- c(arrow_label, only7ds_label, only2021_label)

# --- Plot ---

p <- ggplot() +
  # background
  geom_point(data = bg_df, aes(x = logFC, y = nlog10_fdr),
             colour = "grey90", size = 0.2, alpha = 0.15) +
  # thresholds
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  # arrows
  geom_segment(data = arrows_df,
               aes(x = x0, y = y0, xend = x1, yend = y1, colour = arrow_label),
               arrow = arrow(length = unit(0.8, "mm"), type = "closed"),
               alpha = 0.08, linewidth = 0.2) +
  # 2021 end (open)
  geom_point(data = arrows_df,
             aes(x = x0, y = y0, colour = arrow_label),
             shape = 1, size = 0.6, alpha = 0.15) +
  # 7ds end (filled)
  geom_point(data = arrows_df,
             aes(x = x1, y = y1, colour = arrow_label),
             shape = 16, size = 0.8, alpha = 0.25) +
  # 7ds-only
  geom_point(data = dots_7ds_only,
             aes(x = logFC, y = nlog10_fdr, colour = only7ds_label),
             shape = 17, size = 1.2, alpha = 0.35) +
  # 2021-only
  geom_point(data = dots_2021_only,
             aes(x = logFC, y = nlog10_fdr, colour = only2021_label),
             shape = 15, size = 1.2, alpha = 0.35) +
  # labels
  geom_text_repel(data = top_arrows,
                  aes(x = x1, y = y1, label = symbol),
                  size = 2.5, fontface = "italic", colour = colours[1],
                  max.overlaps = 20, min.segment.length = 0) +
  geom_text_repel(data = top_7ds_only,
                  aes(x = logFC, y = nlog10_fdr, label = symbol),
                  size = 2.5, fontface = "italic", colour = colours[2],
                  max.overlaps = 20, min.segment.length = 0) +
  geom_text_repel(data = top_2021_only,
                  aes(x = logFC, y = nlog10_fdr, label = symbol),
                  size = 2.5, fontface = "italic", colour = colours[3],
                  max.overlaps = 20, min.segment.length = 0) +
  scale_colour_manual(values = colours, name = NULL) +
  labs(x = expression(log[2]~"fold change"),
       y = expression(-log[10]~"FDR"),
       title = "DEG transitions: 2021 (4 Affy, 22 samples) vs 7ds integration (123 samples)",
       subtitle = "DEG = FDR < 0.05 (no logFC filter). Arrows: open circle = 2021, filled = 7ds.") +
  coord_cartesian(xlim = c(-5.5, 4.5), ylim = c(0, fdr_cap + 0.5)) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 11)) +
  guides(colour = guide_legend(ncol = 1,
    override.aes = list(size = 3, alpha = 1, linewidth = 0.8,
                        shape = c(16, 17, 15))))

ggsave(file.path(output_dir, "fig_volcano_arrows_simple_all.pdf"), p, width = 11, height = 9)
png(file.path(output_dir, "fig_volcano_arrows_simple_all.png"),
    width = 3300, height = 2700, res = 300, type = "cairo")
print(p)
invisible(dev.off())

cat("Saved to:", file.path(output_dir, "fig_volcano_arrows_simple_all.{pdf,png}"), "\n")
