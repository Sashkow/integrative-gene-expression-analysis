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
lykhenko_all <- read.csv("output/dissertation/integrative_pipeline_disser/1_2/difexp/difexp_all.csv",
                         stringsAsFactors = FALSE)
lykhenko_sig <- read.csv("output/dissertation/integrative_pipeline_disser/1_2/difexp/difexp_filtered.csv",
                         stringsAsFactors = FALSE)

lykhenko_all$ENTREZID <- as.character(lykhenko_all$ENTREZID)
lykhenko_sig_ids <- as.character(lykhenko_sig$ENTREZID)

gained_ids <- setdiff(as.character(sig_softimpute$gene), as.character(sig_intersection$gene))
intersection_ids <- as.character(sig_intersection$gene)

# --- Build volcano data with categories ---

volcano_df <- data.frame(
  entrezid = lykhenko_all$ENTREZID,
  symbol   = lykhenko_all$SYMBOL,
  logFC    = lykhenko_all$logFC,
  fdr      = lykhenko_all$adj.P.Val,
  stringsAsFactors = FALSE
)
volcano_df$neg_log10_fdr <- -log10(pmax(volcano_df$fdr, 1e-50))

is_gained <- volcano_df$entrezid %in% gained_ids
is_intersection <- volcano_df$entrezid %in% intersection_ids
is_lykhenko_sig <- volcano_df$entrezid %in% lykhenko_sig_ids

volcano_df$group <- "Other"
volcano_df$group[is_intersection & is_lykhenko_sig] <- "Intersection + Lykhenko 2021"
volcano_df$group[is_intersection & !is_lykhenko_sig] <- "Intersection only"
volcano_df$group[is_gained & is_lykhenko_sig] <- "Gained + Lykhenko 2021"
volcano_df$group[is_gained & !is_lykhenko_sig] <- "Gained only"

# --- Counts for legend ---

counts <- table(volcano_df$group)
label_map <- c(
  "Gained + Lykhenko 2021"       = sprintf("Gained + Lykhenko 2021 (n=%d)", counts["Gained + Lykhenko 2021"]),
  "Gained only"                   = sprintf("Gained only (n=%d)", counts["Gained only"]),
  "Intersection + Lykhenko 2021" = sprintf("Intersection + Lykhenko 2021 (n=%d)", counts["Intersection + Lykhenko 2021"]),
  "Intersection only"             = sprintf("Intersection only (n=%d)", counts["Intersection only"]),
  "Other"                         = sprintf("Other (n=%d)", counts["Other"])
)

volcano_df$group_label <- label_map[volcano_df$group]

draw_order <- c("Other", "Intersection only", "Intersection + Lykhenko 2021",
                "Gained only", "Gained + Lykhenko 2021")
volcano_df$group <- factor(volcano_df$group, levels = draw_order)
volcano_df <- volcano_df[order(volcano_df$group), ]

level_labels <- label_map[draw_order]
volcano_df$group_label <- factor(volcano_df$group_label, levels = level_labels)

# --- Colours and sizes ---

colours <- c(
  "Other"                         = "grey80",
  "Intersection only"             = "#92C5DE",
  "Intersection + Lykhenko 2021" = "#2166AC",
  "Gained only"                   = "#FDAE61",
  "Gained + Lykhenko 2021"       = "#D6604D"
)
names(colours) <- level_labels

sizes <- c(0.4, 1.8, 1.8, 2.2, 2.2)
names(sizes) <- level_labels

alphas <- c(0.12, 0.7, 0.8, 0.7, 0.8)
names(alphas) <- level_labels

# --- Label top gained genes ---

top_gained <- volcano_df[volcano_df$group %in% c("Gained + Lykhenko 2021", "Gained only"), ]
top_gained <- top_gained[order(-top_gained$neg_log10_fdr), ][1:12, ]

# --- Plot ---

p <- ggplot(volcano_df, aes(x = logFC, y = neg_log10_fdr, colour = group_label)) +
  geom_point(aes(size = group_label, alpha = group_label)) +
  geom_text_repel(data = top_gained, aes(label = symbol), size = 3,
                  max.overlaps = 20, show.legend = FALSE, fontface = "italic") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = colours, name = NULL) +
  scale_size_manual(values = sizes, name = NULL) +
  scale_alpha_manual(values = alphas, name = NULL) +
  labs(x = expression(log[2]~"fold change (Lykhenko 2021)"),
       y = expression(-log[10]~"FDR (Lykhenko 2021)"),
       title = "Gained DEGs were already detectable in the 2021 analysis") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.margin = margin(t = 0)) +
  guides(colour = guide_legend(ncol = 2, override.aes = list(size = 3, alpha = 1)))

ggsave(file.path(output_dir, "fig_volcano_2021_gained_v2.pdf"), p,
       width = 9, height = 7.5)
png(file.path(output_dir, "fig_volcano_2021_gained_v2.png"),
    width = 2700, height = 2250, res = 300, type = "cairo")
print(p)
invisible(dev.off())

cat("Counts per group:\n")
print(counts)
cat("\nSaved to:", file.path(output_dir, "fig_volcano_2021_gained_v2.{pdf,png}"), "\n")
