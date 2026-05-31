options(bitmapType = "cairo")
library(ggplot2)

output_dir <- "output/article_validation/prior_study_concordance"
base_dir <- "output/yehor_sashko/phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko"

# --- Load full limma tables ---

full_7ds <- read.delim(file.path(base_dir, "difexp_softimpute_combat_ref.tsv"),
                       stringsAsFactors = FALSE)
lykhenko_all <- read.csv("output/dissertation/integrative_pipeline_disser/1_2/difexp/difexp_all.csv",
                         stringsAsFactors = FALSE)

lykhenko_all$ENTREZID <- as.character(lykhenko_all$ENTREZID)
full_7ds$gene <- as.character(full_7ds$gene)

# --- DEG sets: FDR < 0.05 only ---

deg_7ds  <- full_7ds$gene[full_7ds$adj.P.Val < 0.05]
deg_2021 <- lykhenko_all$ENTREZID[lykhenko_all$adj.P.Val < 0.05]

all_degs <- unique(c(deg_7ds, deg_2021))
in_7ds  <- all_degs %in% deg_7ds
in_2021 <- all_degs %in% deg_2021

both_ids      <- all_degs[in_7ds & in_2021]
only_7ds_ids  <- all_degs[in_7ds & !in_2021]
only_2021_ids <- all_degs[!in_7ds & in_2021]

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
  entrezid = both_ids,
  x0 = lykhenko_all$logFC[idx_2021],
  y0 = pmin(-log10(pmax(lykhenko_all$adj.P.Val[idx_2021], 1e-50)), fdr_cap),
  x1 = full_7ds$logFC[idx_7ds],
  y1 = pmin(-log10(pmax(full_7ds$adj.P.Val[idx_7ds], 1e-50)), fdr_cap),
  stringsAsFactors = FALSE
)

# --- 7ds-only and 2021-only dots ---

idx_7ds_only <- match(only_7ds_ids, full_7ds$gene)
dots_7ds_only <- data.frame(
  logFC = full_7ds$logFC[idx_7ds_only],
  nlog10_fdr = pmin(-log10(pmax(full_7ds$adj.P.Val[idx_7ds_only], 1e-50)), fdr_cap),
  stringsAsFactors = FALSE
)

idx_2021_only <- match(only_2021_ids, lykhenko_all$ENTREZID)
dots_2021_only <- data.frame(
  logFC = lykhenko_all$logFC[idx_2021_only],
  nlog10_fdr = pmin(-log10(pmax(lykhenko_all$adj.P.Val[idx_2021_only], 1e-50)), fdr_cap),
  stringsAsFactors = FALSE
)

# --- Cluster arrows ---

feat <- arrows_df[, c("x0", "y0", "x1", "y1")]
feat_scaled <- scale(feat)

best_k <- 10
set.seed(42)
km <- kmeans(feat_scaled, centers = best_k, nstart = 50, iter.max = 100)
arrows_df$cluster <- km$cluster

# Merge tiny clusters (< 20 genes for this larger set)
cl_sizes <- table(arrows_df$cluster)
small_cls <- as.integer(names(cl_sizes[cl_sizes < 20]))

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

# --- Cluster summaries ---

cluster_info <- do.call(rbind, lapply(sort(unique(arrows_df$cluster)), function(cl) {
  sub <- arrows_df[arrows_df$cluster == cl, ]
  pad_x <- 0.06
  pad_y <- 0.12
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

# Colours for non-arrow groups
col_7ds_only  <- "#D6604D"
col_2021_only <- "#5AAE61"

# --- Plot ---

p <- ggplot() +
  # background
  geom_point(data = bg_df, aes(x = logFC, y = nlog10_fdr),
             colour = "grey92", size = 0.15, alpha = 0.15) +
  # thresholds
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey55", linewidth = 0.3) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey55", linewidth = 0.3) +
  # 7ds-only dots
  geom_point(data = dots_7ds_only,
             aes(x = logFC, y = nlog10_fdr),
             colour = col_7ds_only, shape = 17, size = 0.6, alpha = 0.15) +
  # 2021-only dots
  geom_point(data = dots_2021_only,
             aes(x = logFC, y = nlog10_fdr),
             colour = col_2021_only, shape = 15, size = 0.6, alpha = 0.15)

# Rectangles and arrows per cluster
for (i in seq_len(nrow(cluster_info))) {
  ci <- cluster_info[i, ]
  col <- cluster_colours[as.character(ci$cluster)]

  # "from" rectangle (dashed)
  p <- p + annotate("rect",
    xmin = ci$from_xmin, xmax = ci$from_xmax,
    ymin = ci$from_ymin, ymax = ci$from_ymax,
    fill = col, alpha = 0.08, colour = col, linetype = "dashed", linewidth = 0.5)

  # "to" rectangle (solid)
  p <- p + annotate("rect",
    xmin = ci$to_xmin, xmax = ci$to_xmax,
    ymin = ci$to_ymin, ymax = ci$to_ymax,
    fill = col, alpha = 0.12, colour = col, linetype = "solid", linewidth = 0.6)

  # trend arrow
  p <- p + annotate("segment",
    x = ci$from_cx, y = ci$from_cy,
    xend = ci$to_cx, yend = ci$to_cy,
    colour = col, linewidth = 1.4, alpha = 0.85,
    arrow = arrow(length = unit(3.5, "mm"), type = "closed"))

  # gene count label at midpoint
  mid_x <- (ci$from_cx + ci$to_cx) / 2
  mid_y <- (ci$from_cy + ci$to_cy) / 2
  # offset label slightly to avoid overlap with arrow
  dx <- ci$to_cx - ci$from_cx
  dy <- ci$to_cy - ci$from_cy
  len <- sqrt(dx^2 + dy^2)
  if (len > 0) {
    off_x <- -dy / len * 0.35
    off_y <-  dx / len * 0.35
  } else {
    off_x <- 0.3; off_y <- 0
  }
  p <- p + annotate("label",
    x = mid_x + off_x, y = mid_y + off_y,
    label = sprintf("n=%d", ci$n),
    size = 3, fontface = "bold", colour = col,
    fill = "white", alpha = 0.9,
    label.padding = unit(1.5, "pt"))
}

# Legend annotation (manual, since we have mixed annotation + geom)
p <- p +
  annotate("text", x = -5.3, y = fdr_cap - 0.3,
    label = sprintf("Clustered arrows: DEG in both (total n=%d)", length(both_ids)),
    hjust = 0, size = 3.2, fontface = "bold") +
  annotate("point", x = -5.3, y = fdr_cap - 1.0,
    colour = col_7ds_only, shape = 17, size = 2.5) +
  annotate("text", x = -5.0, y = fdr_cap - 1.0,
    label = sprintf("DEG only in 7ds (n=%d)", length(only_7ds_ids)),
    hjust = 0, size = 3, colour = col_7ds_only) +
  annotate("point", x = -5.3, y = fdr_cap - 1.7,
    colour = col_2021_only, shape = 15, size = 2.5) +
  annotate("text", x = -5.0, y = fdr_cap - 1.7,
    label = sprintf("DEG only in 2021 (n=%d)", length(only_2021_ids)),
    hjust = 0, size = 3, colour = col_2021_only)

p <- p +
  labs(x = expression(log[2]~"fold change"),
       y = expression(-log[10]~"FDR"),
       title = "Clustered DEG transitions: 2021 (4 Affy, 22 samples) -> 7ds (123 samples)",
       subtitle = "DEG = FDR < 0.05. Dashed rect = 2021 region, solid rect = 7ds region, arrows = cluster trend.") +
  coord_cartesian(xlim = c(-5.5, 4.5), ylim = c(0, fdr_cap + 0.5)) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        plot.title = element_text(size = 11))

ggsave(file.path(output_dir, "fig_volcano_clustered_all.pdf"), p, width = 11, height = 9)
png(file.path(output_dir, "fig_volcano_clustered_all.png"),
    width = 3300, height = 2700, res = 300, type = "cairo")
print(p)
invisible(dev.off())

cat("\nCluster details:\n")
for (i in seq_len(nrow(cluster_info))) {
  ci <- cluster_info[i, ]
  cat(sprintf("  Cluster %2d (n=%4d): from (%5.1f,%4.1f) -> to (%5.1f,%4.1f)  delta=(%5.1f,%5.1f)\n",
    ci$cluster, ci$n, ci$from_cx, ci$from_cy, ci$to_cx, ci$to_cy,
    ci$to_cx - ci$from_cx, ci$to_cy - ci$from_cy))
}

cat("\nSaved to:", file.path(output_dir, "fig_volcano_clustered_all.{pdf,png}"), "\n")
