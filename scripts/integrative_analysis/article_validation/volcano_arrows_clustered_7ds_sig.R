options(bitmapType = "cairo")
library(ggplot2)

output_dir <- "output/article_validation/prior_study_concordance"
base_dir <- "output/yehor_sashko/phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko"

# --- Load data ---

full_7ds <- read.delim(file.path(base_dir, "difexp_softimpute_combat_ref.tsv"),
                       stringsAsFactors = FALSE)
lykhenko_all <- read.csv("output/dissertation/integrative_pipeline_disser/1_2/difexp/difexp_all.csv",
                         stringsAsFactors = FALSE)
lykhenko_sig <- read.csv("output/dissertation/integrative_pipeline_disser/1_2/difexp/difexp_filtered.csv",
                         stringsAsFactors = FALSE)

lykhenko_all$ENTREZID <- as.character(lykhenko_all$ENTREZID)
lykhenko_sig$ENTREZID <- as.character(lykhenko_sig$ENTREZID)
full_7ds$gene <- as.character(full_7ds$gene)

fdr_cap <- 14

# --- Arrows: all 7ds FDR<0.05 genes that exist in 2021 table ---

deg_7ds_ids <- full_7ds$gene[full_7ds$adj.P.Val < 0.05]
arrow_ids <- deg_7ds_ids[deg_7ds_ids %in% lykhenko_all$ENTREZID]
no_2021_ids <- deg_7ds_ids[!deg_7ds_ids %in% lykhenko_all$ENTREZID]

# 2021-only: FDR<0.05 in 2021 but NOT in 7ds
deg_2021_ids <- lykhenko_all$ENTREZID[lykhenko_all$adj.P.Val < 0.05]
only_2021_ids <- setdiff(deg_2021_ids, deg_7ds_ids)
only_2021_ids <- only_2021_ids[only_2021_ids %in% lykhenko_all$ENTREZID]

cat("Arrows (7ds FDR<0.05, in 2021 table):", length(arrow_ids), "\n")
cat("7ds-only (no 2021 data):", length(no_2021_ids), "\n")
cat("2021-only (FDR<0.05 in 2021, not in 7ds):", length(only_2021_ids), "\n\n")

# --- Build arrow data ---

idx_2021 <- match(arrow_ids, lykhenko_all$ENTREZID)
idx_7ds  <- match(arrow_ids, full_7ds$gene)

arrows_df <- data.frame(
  entrezid = arrow_ids,
  x0 = lykhenko_all$logFC[idx_2021],
  y0 = pmin(-log10(pmax(lykhenko_all$adj.P.Val[idx_2021], 1e-50)), fdr_cap),
  x1 = full_7ds$logFC[idx_7ds],
  y1 = pmin(-log10(pmax(full_7ds$adj.P.Val[idx_7ds], 1e-50)), fdr_cap),
  stringsAsFactors = FALSE
)

# --- 7ds-only dots (no 2021 position) ---

idx_no2021 <- match(no_2021_ids, full_7ds$gene)
dots_no2021 <- data.frame(
  logFC = full_7ds$logFC[idx_no2021],
  nlog10_fdr = pmin(-log10(pmax(full_7ds$adj.P.Val[idx_no2021], 1e-50)), fdr_cap),
  stringsAsFactors = FALSE
)

# --- 2021-only dots ---

idx_only2021 <- match(only_2021_ids, lykhenko_all$ENTREZID)
dots_only2021 <- data.frame(
  logFC = lykhenko_all$logFC[idx_only2021],
  nlog10_fdr = pmin(-log10(pmax(lykhenko_all$adj.P.Val[idx_only2021], 1e-50)), fdr_cap),
  stringsAsFactors = FALSE
)

# --- Cluster arrows ---

feat <- arrows_df[, c("x0", "y0", "x1", "y1")]
feat_scaled <- scale(feat)

set.seed(42)
km <- kmeans(feat_scaled, centers = 6, nstart = 50, iter.max = 100)
arrows_df$cluster <- km$cluster

# Merge small clusters (< 30 genes)
cl_sizes <- table(arrows_df$cluster)
small_cls <- as.integer(names(cl_sizes[cl_sizes < 30]))

if (length(small_cls) > 0) {
  centers <- km$centers
  for (sc in small_cls) {
    remaining <- setdiff(seq_len(nrow(centers)), small_cls)
    if (length(remaining) == 0) break
    dists <- apply(centers[remaining, , drop = FALSE], 1, function(c)
      sum((centers[sc, ] - c)^2))
    nearest <- remaining[which.min(dists)]
    arrows_df$cluster[arrows_df$cluster == sc] <- nearest
  }
  arrows_df$cluster <- as.integer(factor(arrows_df$cluster))
}

n_clusters <- length(unique(arrows_df$cluster))
cat("Final clusters:", n_clusters, "\n")
cat("Cluster sizes:\n")
print(table(arrows_df$cluster))

# --- Cluster summaries ---

# Use 2.5th and 97.5th percentiles for rectangle bounds (trim 5% outliers)
cluster_info <- do.call(rbind, lapply(sort(unique(arrows_df$cluster)), function(cl) {
  sub <- arrows_df[arrows_df$cluster == cl, ]
  pad_x <- 0.05
  pad_y <- 0.1
  data.frame(
    cluster = cl,
    n = nrow(sub),
    from_xmin = quantile(sub$x0, 0.025) - pad_x,
    from_xmax = quantile(sub$x0, 0.975) + pad_x,
    from_ymin = quantile(sub$y0, 0.025) - pad_y,
    from_ymax = quantile(sub$y0, 0.975) + pad_y,
    from_cx   = mean(sub$x0),
    from_cy   = mean(sub$y0),
    to_xmin   = quantile(sub$x1, 0.025) - pad_x,
    to_xmax   = quantile(sub$x1, 0.975) + pad_x,
    to_ymin   = quantile(sub$y1, 0.025) - pad_y,
    to_ymax   = quantile(sub$y1, 0.975) + pad_y,
    to_cx     = mean(sub$x1),
    to_cy     = mean(sub$y1),
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

col_no2021   <- "#D6604D"
col_only2021 <- "#5AAE61"

# --- Plot ---

p <- ggplot() +
  geom_point(data = bg_df, aes(x = logFC, y = nlog10_fdr),
             colour = "grey92", size = 0.15, alpha = 0.12) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey55", linewidth = 0.3) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey55", linewidth = 0.3) +
  # 7ds-only dots (no 2021 data)
  geom_point(data = dots_no2021,
             aes(x = logFC, y = nlog10_fdr),
             colour = col_no2021, shape = 17, size = 1.5, alpha = 0.5) +
  # 2021-only dots
  geom_point(data = dots_only2021,
             aes(x = logFC, y = nlog10_fdr),
             colour = col_only2021, shape = 15, size = 0.5, alpha = 0.12)

# Rectangles and arrows
for (i in seq_len(nrow(cluster_info))) {
  ci <- cluster_info[i, ]
  col <- cluster_colours[as.character(ci$cluster)]

  p <- p + annotate("rect",
    xmin = ci$from_xmin, xmax = ci$from_xmax,
    ymin = ci$from_ymin, ymax = ci$from_ymax,
    fill = col, alpha = 0.07, colour = col, linetype = "dashed", linewidth = 0.5)

  p <- p + annotate("rect",
    xmin = ci$to_xmin, xmax = ci$to_xmax,
    ymin = ci$to_ymin, ymax = ci$to_ymax,
    fill = col, alpha = 0.11, colour = col, linetype = "solid", linewidth = 0.6)

  p <- p + annotate("segment",
    x = ci$from_cx, y = ci$from_cy,
    xend = ci$to_cx, yend = ci$to_cy,
    colour = col, linewidth = 1.5, alpha = 0.85,
    arrow = arrow(length = unit(3.5, "mm"), type = "closed"))

  # label offset perpendicular to arrow
  dx <- ci$to_cx - ci$from_cx
  dy <- ci$to_cy - ci$from_cy
  len <- sqrt(dx^2 + dy^2)
  if (len > 0.3) {
    off_x <- -dy / len * 0.4
    off_y <-  dx / len * 0.4
  } else {
    off_x <- 0.35; off_y <- 0
  }
  mid_x <- (ci$from_cx + ci$to_cx) / 2 + off_x
  mid_y <- (ci$from_cy + ci$to_cy) / 2 + off_y

  p <- p + annotate("label",
    x = mid_x, y = mid_y,
    label = sprintf("n=%d", ci$n),
    size = 3, fontface = "bold", colour = col,
    fill = "white", alpha = 0.9,
    label.padding = unit(1.5, "pt"))
}

# --- Build cluster legend descriptions ---

cluster_info$direction <- ifelse(cluster_info$to_cx - cluster_info$from_cx > 0.2, "up-regulated shift",
                          ifelse(cluster_info$to_cx - cluster_info$from_cx < -0.2, "down-regulated shift",
                          "stable logFC"))
cluster_info$strength <- ifelse(cluster_info$to_cy - cluster_info$from_cy > 5, "strong FDR gain",
                         ifelse(cluster_info$to_cy - cluster_info$from_cy > 2, "moderate FDR gain",
                         "modest FDR gain"))

# Build legend data frame for geom_point + geom_segment in aes so ggplot makes a real legend
legend_df <- data.frame(
  cluster_label = sprintf("Cluster %d: n=%d, %s, %s",
                          cluster_info$cluster, cluster_info$n,
                          cluster_info$direction, cluster_info$strength),
  col = cluster_colours[as.character(cluster_info$cluster)],
  stringsAsFactors = FALSE
)

# Add dummy geom layers that produce legend entries for each cluster
for (i in seq_len(nrow(cluster_info))) {
  ci <- cluster_info[i, ]
  col <- cluster_colours[as.character(ci$cluster)]
  lbl <- legend_df$cluster_label[i]

  p <- p + geom_segment(
    data = data.frame(x = ci$from_cx, y = ci$from_cy,
                      xend = ci$to_cx, yend = ci$to_cy,
                      label = lbl, stringsAsFactors = FALSE),
    aes(x = x, y = y, xend = xend, yend = yend, colour = label),
    arrow = arrow(length = unit(2.5, "mm"), type = "closed"),
    linewidth = 1.2, alpha = 0, show.legend = TRUE)
}

# Add non-arrow legend entries as dummy points
p <- p +
  geom_point(data = data.frame(x = -99, y = -99,
    label = sprintf("7ds FDR<0.05, not in 2021 table (n=%d)", length(no_2021_ids)),
    stringsAsFactors = FALSE),
    aes(x = x, y = y, colour = label), shape = 17, size = 2.5, show.legend = TRUE) +
  geom_point(data = data.frame(x = -99, y = -99,
    label = sprintf("2021 FDR<0.05, not in 7ds (n=%d)", length(only_2021_ids)),
    stringsAsFactors = FALSE),
    aes(x = x, y = y, colour = label), shape = 15, size = 2.5, show.legend = TRUE)

# Colour scale with all entries
all_legend_colours <- c(
  setNames(cluster_colours[as.character(cluster_info$cluster)], legend_df$cluster_label),
  setNames(col_no2021, sprintf("7ds FDR<0.05, not in 2021 table (n=%d)", length(no_2021_ids))),
  setNames(col_only2021, sprintf("2021 FDR<0.05, not in 7ds (n=%d)", length(only_2021_ids)))
)

p <- p +
  scale_colour_manual(values = all_legend_colours, name = NULL) +
  labs(x = expression(log[2]~"fold change"),
       y = expression(-log[10]~"FDR"),
       title = "All 7ds-significant genes: clustered transitions from 2021 to 7-dataset integration",
       subtitle = "Arrows for every gene with FDR<0.05 in 7ds, regardless of 2021 significance.\nDashed rect = 2021 region (95% of cluster), solid rect = 7ds region.") +
  coord_cartesian(xlim = c(-5.5, 4.5), ylim = c(0, fdr_cap + 0.5)) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 8),
        legend.key.width = unit(15, "pt"),
        plot.title = element_text(size = 11)) +
  guides(colour = guide_legend(ncol = 2, override.aes = list(alpha = 1, linewidth = 1, size = 2.5)))

ggsave(file.path(output_dir, "fig_volcano_clustered_7ds_sig.pdf"), p, width = 11, height = 9)
png(file.path(output_dir, "fig_volcano_clustered_7ds_sig.png"),
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
cat("\nSaved to:", file.path(output_dir, "fig_volcano_clustered_7ds_sig.{pdf,png}"), "\n")
