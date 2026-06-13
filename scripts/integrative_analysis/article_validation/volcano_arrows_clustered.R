options(bitmapType = "cairo")
library(ggplot2)

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

all_7ds_deg_ids <- unique(c(as.character(sig_softimpute$gene),
                            as.character(sig_intersection$gene)))
lykhenko_sig_ids <- lykhenko_sig$ENTREZID

both_ids <- intersect(all_7ds_deg_ids, lykhenko_sig_ids)
both_ids <- both_ids[both_ids %in% lykhenko_all$ENTREZID & both_ids %in% full_7ds$gene]

idx_2021 <- match(both_ids, lykhenko_all$ENTREZID)
idx_7ds  <- match(both_ids, full_7ds$gene)

fdr_cap <- 14

arrows_df <- data.frame(
  entrezid   = both_ids,
  x0 = lykhenko_all$logFC[idx_2021],
  y0 = pmin(-log10(pmax(lykhenko_all$adj.P.Val[idx_2021], 1e-50)), fdr_cap),
  x1 = full_7ds$logFC[idx_7ds],
  y1 = pmin(-log10(pmax(full_7ds$adj.P.Val[idx_7ds], 1e-50)), fdr_cap),
  stringsAsFactors = FALSE
)

# --- Cluster on (from_x, from_y, to_x, to_y) ---

# Normalize each dimension to [0,1] for balanced clustering
feat <- arrows_df[, c("x0", "y0", "x1", "y1")]
feat_scaled <- scale(feat)

# Fixed k — 8 gives a good balance of readability and detail
best_k <- 8

set.seed(42)
km <- kmeans(feat_scaled, centers = best_k, nstart = 50, iter.max = 100)
arrows_df$cluster <- km$cluster

# --- Remove tiny clusters (< 3 genes): merge into nearest ---

cl_sizes <- table(arrows_df$cluster)
small_cls <- as.integer(names(cl_sizes[cl_sizes < 3]))

if (length(small_cls) > 0) {
  centers <- km$centers
  for (sc in small_cls) {
    dists <- apply(centers[-sc, , drop = FALSE], 1, function(c)
      sum((centers[sc, ] - c)^2))
    nearest <- as.integer(names(which.min(dists)))
    arrows_df$cluster[arrows_df$cluster == sc] <- nearest
  }
  arrows_df$cluster <- as.integer(factor(arrows_df$cluster))
}

n_clusters <- length(unique(arrows_df$cluster))
cat("Final clusters:", n_clusters, "\n")
cat("Cluster sizes:\n")
print(table(arrows_df$cluster))

# --- Compute cluster summaries ---

cluster_info <- do.call(rbind, lapply(sort(unique(arrows_df$cluster)), function(cl) {
  sub <- arrows_df[arrows_df$cluster == cl, ]
  pad_x <- 0.08
  pad_y <- 0.15
  data.frame(
    cluster = cl,
    n = nrow(sub),
    from_xmin = min(sub$x0) - pad_x, from_xmax = max(sub$x0) + pad_x,
    from_ymin = min(sub$y0) - pad_y, from_ymax = max(sub$y0) + pad_y,
    from_cx   = mean(sub$x0),        from_cy   = mean(sub$y0),
    to_xmin   = min(sub$x1) - pad_x, to_xmax   = max(sub$x1) + pad_x,
    to_ymin   = min(sub$y1) - pad_y, to_ymax   = max(sub$y1) + pad_y,
    to_cx     = mean(sub$x1),        to_cy     = mean(sub$y1),
    stringsAsFactors = FALSE
  )
}))

# --- Background ---

bg_df <- data.frame(
  logFC = lykhenko_all$logFC,
  nlog10_fdr = pmin(-log10(pmax(lykhenko_all$adj.P.Val, 1e-50)), fdr_cap),
  stringsAsFactors = FALSE
)

# --- Palette ---

pal <- c("#2166AC", "#D6604D", "#5AAE61", "#F4A582", "#8073AC",
         "#FDB863", "#1B9E77", "#E7298A", "#66A61E", "#7570B3",
         "#E6AB02", "#A6761D", "#E41A1C", "#377EB8", "#984EA3")
cluster_colours <- pal[seq_len(n_clusters)]
names(cluster_colours) <- as.character(sort(unique(arrows_df$cluster)))

# --- Legend labels ---

cluster_info$label <- sprintf("Cluster %d (n=%d)", cluster_info$cluster, cluster_info$n)

# --- Plot ---

p <- ggplot() +
  # background
  geom_point(data = bg_df, aes(x = logFC, y = nlog10_fdr),
             colour = "grey90", size = 0.2, alpha = 0.2) +
  # threshold lines
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey60", linewidth = 0.3)

# Add "from" and "to" rectangles + arrows per cluster
for (i in seq_len(nrow(cluster_info))) {
  ci <- cluster_info[i, ]
  col <- cluster_colours[as.character(ci$cluster)]

  # "from" rectangle (dashed border)
  p <- p + annotate("rect",
    xmin = ci$from_xmin, xmax = ci$from_xmax,
    ymin = ci$from_ymin, ymax = ci$from_ymax,
    fill = col, alpha = 0.10, colour = col, linetype = "dashed", linewidth = 0.5)

  # "to" rectangle (solid border)
  p <- p + annotate("rect",
    xmin = ci$to_xmin, xmax = ci$to_xmax,
    ymin = ci$to_ymin, ymax = ci$to_ymax,
    fill = col, alpha = 0.15, colour = col, linetype = "solid", linewidth = 0.6)

  # arrow from center of "from" to center of "to"
  p <- p + annotate("segment",
    x = ci$from_cx, y = ci$from_cy,
    xend = ci$to_cx, yend = ci$to_cy,
    colour = col, linewidth = 1.2, alpha = 0.8,
    arrow = arrow(length = unit(3, "mm"), type = "closed"))

  # label: gene count at the midpoint of the arrow
  mid_x <- (ci$from_cx + ci$to_cx) / 2
  mid_y <- (ci$from_cy + ci$to_cy) / 2
  p <- p + annotate("label",
    x = mid_x, y = mid_y,
    label = sprintf("n=%d", ci$n),
    size = 3, fontface = "bold", colour = col,
    fill = "white", alpha = 0.85, label.size = 0.3,
    label.padding = unit(1.5, "pt"))
}

p <- p +
  labs(x = expression(log[2]~"fold change"),
       y = expression(-log[10]~"FDR"),
       title = "Clustered DEG transitions: 2021 (4 Affy, 22 samples) → 7ds (123 samples)",
       subtitle = "Dashed rectangles = 2021 region, solid rectangles = 7ds region, arrows = trend direction") +
  coord_cartesian(xlim = c(-5.5, 4.5), ylim = c(0, fdr_cap + 0.5)) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        plot.title = element_text(size = 11))

ggsave(file.path(output_dir, "fig_volcano_arrows_clustered.pdf"), p, width = 11, height = 9)
png(file.path(output_dir, "fig_volcano_arrows_clustered.png"),
    width = 3300, height = 2700, res = 300, type = "cairo")
print(p)
invisible(dev.off())

cat("\nCluster details:\n")
for (i in seq_len(nrow(cluster_info))) {
  ci <- cluster_info[i, ]
  cat(sprintf("  Cluster %d (n=%2d): from (%.1f,%.1f) → to (%.1f,%.1f)  Δ=(%.1f,%.1f)\n",
    ci$cluster, ci$n, ci$from_cx, ci$from_cy, ci$to_cx, ci$to_cy,
    ci$to_cx - ci$from_cx, ci$to_cy - ci$from_cy))
}

cat("\nSaved to:", file.path(output_dir, "fig_volcano_arrows_clustered.{pdf,png}"), "\n")
