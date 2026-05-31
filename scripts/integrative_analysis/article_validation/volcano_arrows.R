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
full_7ds_intersect <- read.delim(file.path(base_dir, "difexp_none_combat_ref.tsv"),
                                 stringsAsFactors = FALSE)
lykhenko_all <- read.csv("output/dissertation/integrative_pipeline_disser/1_2/difexp/difexp_all.csv",
                         stringsAsFactors = FALSE)

lykhenko_all$ENTREZID <- as.character(lykhenko_all$ENTREZID)
full_7ds$gene <- as.character(full_7ds$gene)
full_7ds_intersect$gene <- as.character(full_7ds_intersect$gene)

gained_ids <- setdiff(as.character(sig_softimpute$gene), as.character(sig_intersection$gene))
intersection_ids <- as.character(sig_intersection$gene)

# --- Build arrow data for GAINED DEGs ---

gained_matched <- gained_ids[gained_ids %in% lykhenko_all$ENTREZID]
idx_2021 <- match(gained_matched, lykhenko_all$ENTREZID)
idx_7ds  <- match(gained_matched, full_7ds$gene)

gained_arrows <- data.frame(
  entrezid   = gained_matched,
  symbol     = lykhenko_all$SYMBOL[idx_2021],
  logFC_2021 = lykhenko_all$logFC[idx_2021],
  fdr_2021   = lykhenko_all$adj.P.Val[idx_2021],
  logFC_7ds  = full_7ds$logFC[idx_7ds],
  fdr_7ds    = full_7ds$adj.P.Val[idx_7ds],
  deg_set    = "Gained",
  stringsAsFactors = FALSE
)

# --- Build arrow data for INTERSECTION DEGs ---

intersect_matched <- intersection_ids[intersection_ids %in% lykhenko_all$ENTREZID]
idx_2021_i <- match(intersect_matched, lykhenko_all$ENTREZID)
idx_7ds_i  <- match(intersect_matched, full_7ds_intersect$gene)

intersect_arrows <- data.frame(
  entrezid   = intersect_matched,
  symbol     = lykhenko_all$SYMBOL[idx_2021_i],
  logFC_2021 = lykhenko_all$logFC[idx_2021_i],
  fdr_2021   = lykhenko_all$adj.P.Val[idx_2021_i],
  logFC_7ds  = full_7ds_intersect$logFC[idx_7ds_i],
  fdr_7ds    = full_7ds_intersect$adj.P.Val[idx_7ds_i],
  deg_set    = "Intersection",
  stringsAsFactors = FALSE
)

arrows_df <- rbind(gained_arrows, intersect_arrows)

arrows_df$nlog10_fdr_2021 <- -log10(pmax(arrows_df$fdr_2021, 1e-50))
arrows_df$nlog10_fdr_7ds  <- -log10(pmax(arrows_df$fdr_7ds, 1e-50))

# --- Categorize transitions ---

arrows_df$was_sig_2021 <- arrows_df$fdr_2021 < 0.05 & abs(arrows_df$logFC_2021) > 1

arrows_df$transition <- ifelse(
  arrows_df$deg_set == "Intersection" & arrows_df$was_sig_2021,
    "Intersection, already sig in 2021",
  ifelse(arrows_df$deg_set == "Intersection",
    "Intersection, newly sig in 7ds",
  ifelse(arrows_df$was_sig_2021,
    "Gained, already sig in 2021",
  ifelse(arrows_df$fdr_2021 < 0.05,
    "Gained, had FDR<0.05 in 2021",
    "Gained, sub-threshold in 2021"))))

cat("Transition counts:\n")
print(table(arrows_df$transition))

# --- Background: all 2021 genes ---

bg_df <- data.frame(
  logFC = lykhenko_all$logFC,
  nlog10_fdr = -log10(pmax(lykhenko_all$adj.P.Val, 1e-50)),
  stringsAsFactors = FALSE
)

# --- Cap -log10 FDR ---

fdr_cap <- 14
bg_df$nlog10_fdr <- pmin(bg_df$nlog10_fdr, fdr_cap)
arrows_df$nlog10_fdr_2021 <- pmin(arrows_df$nlog10_fdr_2021, fdr_cap)
arrows_df$nlog10_fdr_7ds  <- pmin(arrows_df$nlog10_fdr_7ds, fdr_cap)

# --- Legend labels with counts ---

trans_counts <- table(arrows_df$transition)
trans_label_map <- sprintf("%s (n=%d)", names(trans_counts), trans_counts)
names(trans_label_map) <- names(trans_counts)

arrows_df$trans_label <- trans_label_map[arrows_df$transition]

draw_order <- c("Intersection, already sig in 2021",
                "Intersection, newly sig in 7ds",
                "Gained, already sig in 2021",
                "Gained, had FDR<0.05 in 2021",
                "Gained, sub-threshold in 2021")
label_order <- trans_label_map[draw_order]
arrows_df$trans_label <- factor(arrows_df$trans_label, levels = label_order)

# --- Draw order: intersection first (behind), gained on top ---

arrows_df <- arrows_df[order(arrows_df$trans_label), ]

# --- Colours ---

trans_colours <- c(
  "#2166AC",   # Intersection, already sig 2021
  "#92C5DE",   # Intersection, newly sig
  "#B2182B",   # Gained, already sig 2021
  "#F4A582",   # Gained, had FDR<0.05
  "#FDAE61"    # Gained, sub-threshold
)
names(trans_colours) <- label_order

# --- Top movers to label ---

arrows_df$delta <- sqrt((arrows_df$logFC_7ds - arrows_df$logFC_2021)^2 +
                        (arrows_df$nlog10_fdr_7ds - arrows_df$nlog10_fdr_2021)^2)

top_gained <- arrows_df[arrows_df$deg_set == "Gained", ]
top_gained <- top_gained[order(-top_gained$delta), ][1:10, ]

top_intersect <- arrows_df[arrows_df$deg_set == "Intersection", ]
top_intersect <- top_intersect[order(-top_intersect$delta), ][1:5, ]

top_movers <- rbind(top_gained, top_intersect)

# --- Plot ---

p <- ggplot() +
  geom_point(data = bg_df, aes(x = logFC, y = nlog10_fdr),
             colour = "grey88", size = 0.3, alpha = 0.25) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey50") +
  # arrows
  geom_segment(data = arrows_df,
               aes(x = logFC_2021, y = nlog10_fdr_2021,
                   xend = logFC_7ds, yend = nlog10_fdr_7ds,
                   colour = trans_label),
               arrow = arrow(length = unit(1.2, "mm"), type = "closed"),
               alpha = 0.30, linewidth = 0.35) +
  # 2021 positions (open circles)
  geom_point(data = arrows_df,
             aes(x = logFC_2021, y = nlog10_fdr_2021, colour = trans_label),
             shape = 1, size = 1.3, alpha = 0.45) +
  # 7ds positions (filled circles)
  geom_point(data = arrows_df,
             aes(x = logFC_7ds, y = nlog10_fdr_7ds, colour = trans_label),
             shape = 16, size = 1.8, alpha = 0.65) +
  # labels
  geom_text_repel(data = top_movers,
                  aes(x = logFC_7ds, y = nlog10_fdr_7ds, label = symbol),
                  size = 2.8, fontface = "italic",
                  max.overlaps = 25, min.segment.length = 0) +
  scale_colour_manual(values = trans_colours, name = NULL) +
  labs(x = expression(log[2]~"fold change"),
       y = expression(-log[10]~"FDR"),
       title = "DEG transitions from 2021 (4 Affy, 22 samples) to 7-dataset integration (123 samples)",
       subtitle = expression("Open circles = 2021 position, filled circles = 7ds position,"~
                             "arrows = transition")) +
  coord_cartesian(xlim = c(-5.5, 4.5), ylim = c(0, fdr_cap + 0.5)) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 9),
        plot.title = element_text(size = 11)) +
  guides(colour = guide_legend(ncol = 2, override.aes = list(size = 3, alpha = 1, linewidth = 0.8)))

ggsave(file.path(output_dir, "fig_volcano_arrows_v2.pdf"), p, width = 11, height = 9)
png(file.path(output_dir, "fig_volcano_arrows_v2.png"),
    width = 3300, height = 2700, res = 300, type = "cairo")
print(p)
invisible(dev.off())

cat("\nSaved to:", file.path(output_dir, "fig_volcano_arrows_v2.{pdf,png}"), "\n")
