options(bitmapType = "cairo")
library(ggplot2)
library(ggrepel)

output_dir <- "output/article_validation/internal_concordance"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

base_dir <- "output/yehor_sashko/phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko"

# --- Load full limma tables ---

full_softimpute <- read.delim(file.path(base_dir, "difexp_softimpute_combat_ref.tsv"),
                              stringsAsFactors = FALSE)
full_intersect  <- read.delim(file.path(base_dir, "difexp_none_combat_ref.tsv"),
                              stringsAsFactors = FALSE)
sig_softimpute  <- read.delim(file.path(base_dir, "difexp_significant_softimpute_combat_ref.tsv"),
                              stringsAsFactors = FALSE)
sig_intersect   <- read.delim(file.path(base_dir, "difexp_significant_none_combat_ref.tsv"),
                              stringsAsFactors = FALSE)

full_softimpute$gene <- as.character(full_softimpute$gene)
full_intersect$gene  <- as.character(full_intersect$gene)

shared_genes <- intersect(full_softimpute$gene, full_intersect$gene)
cat("Shared genes (intersection set):", length(shared_genes), "\n")
cat("softImpute-only genes:", sum(!full_softimpute$gene %in% full_intersect$gene), "\n")
cat("softImpute DEGs:", nrow(sig_softimpute), "\n")
cat("Intersection DEGs:", nrow(sig_intersect), "\n\n")

# --- Build comparison table for shared genes ---

idx_si <- match(shared_genes, full_softimpute$gene)
idx_io <- match(shared_genes, full_intersect$gene)

comp <- data.frame(
  gene           = shared_genes,
  logFC_softimpute = full_softimpute$logFC[idx_si],
  fdr_softimpute   = full_softimpute$adj.P.Val[idx_si],
  logFC_intersect  = full_intersect$logFC[idx_io],
  fdr_intersect    = full_intersect$adj.P.Val[idx_io],
  stringsAsFactors = FALSE
)

comp$same_direction <- sign(comp$logFC_softimpute) == sign(comp$logFC_intersect)

classify <- function(fdr, logfc) {
  ifelse(fdr < 0.05 & abs(logfc) > 1, "DEG",
  ifelse(fdr < 0.05, "FDR-only",
  ifelse(abs(logfc) > 1, "logFC-only", "NS")))
}

comp$state_softimpute <- classify(comp$fdr_softimpute, comp$logFC_softimpute)
comp$state_intersect  <- classify(comp$fdr_intersect, comp$logFC_intersect)

r_val <- cor(comp$logFC_softimpute, comp$logFC_intersect, use = "complete.obs")
cat(sprintf("logFC Pearson r (shared genes): %.4f\n", r_val))
cat(sprintf("Same direction: %d/%d (%.1f%%)\n",
    sum(comp$same_direction), nrow(comp),
    100 * sum(comp$same_direction) / nrow(comp)))

# --- Transition matrix ---

lvls <- c("DEG", "FDR-only", "logFC-only", "NS")
comp$state_intersect  <- factor(comp$state_intersect, levels = lvls)
comp$state_softimpute <- factor(comp$state_softimpute, levels = lvls)

tab <- table("Intersection_only" = comp$state_intersect,
             "SoftImpute" = comp$state_softimpute)

cat("\nTransition matrix (shared 8,260 genes):\n")
print(tab)
write.csv(as.data.frame.matrix(tab), file.path(output_dir, "transition_matrix_internal.csv"))

# --- Transition matrix heatmap ---

heat_df <- as.data.frame(tab)
names(heat_df) <- c("Origin", "Destination", "Count")
heat_df$Origin <- factor(heat_df$Origin, levels = rev(lvls))
heat_df$Destination <- factor(heat_df$Destination, levels = lvls)
heat_df$pct <- ave(heat_df$Count, heat_df$Origin, FUN = function(x) x / sum(x) * 100)
heat_df$label <- ifelse(heat_df$Count > 0,
  sprintf("%s\n(%.0f%%)", formatC(heat_df$Count, big.mark = ","), heat_df$pct), "")

p_heat <- ggplot(heat_df, aes(x = Destination, y = Origin, fill = log10(Count + 1))) +
  geom_tile(colour = "white", linewidth = 1.2) +
  geom_text(aes(label = label), size = 3.5, lineheight = 0.85) +
  scale_fill_gradient2(low = "white", mid = "#FDDBC7", high = "#2166AC",
                       midpoint = 2, na.value = "grey95",
                       name = expression(log[10]~"(count+1)")) +
  labs(x = "SoftImpute + ComBat-ref (17,531 genes tested)",
       y = "Intersection-only + ComBat-ref (8,260 genes tested)",
       title = "Gene state transitions: intersection-only vs softImpute (shared 8,260 genes)",
       subtitle = "DEG = FDR<0.05 & |logFC|>1, FDR-only = FDR<0.05 & |logFC|≤1, NS = not significant") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 11),
        axis.text.y = element_text(size = 11),
        panel.grid = element_blank(),
        legend.position = "right",
        plot.title = element_text(size = 12, face = "bold"))

ggsave(file.path(output_dir, "fig_transition_matrix_internal.pdf"), p_heat, width = 9, height = 7)
png(file.path(output_dir, "fig_transition_matrix_internal.png"),
    width = 2700, height = 2100, res = 300, type = "cairo")
print(p_heat)
invisible(dev.off())
cat("Saved transition matrix heatmap\n")

# --- logFC scatter plot (shared genes, DEGs highlighted) ---

comp$deg_group <- ifelse(
  comp$gene %in% sig_softimpute$gene & comp$gene %in% sig_intersect$gene, "DEG in both",
  ifelse(comp$gene %in% sig_softimpute$gene, "DEG in softImpute only",
  ifelse(comp$gene %in% sig_intersect$gene, "DEG in intersection only", "Not DEG")))

comp$deg_group <- factor(comp$deg_group,
  levels = c("DEG in both", "DEG in softImpute only", "DEG in intersection only", "Not DEG"))

p_scatter <- ggplot(comp, aes(x = logFC_intersect, y = logFC_softimpute, colour = deg_group)) +
  geom_point(data = comp[comp$deg_group == "Not DEG", ], alpha = 0.05, size = 0.5) +
  geom_point(data = comp[comp$deg_group != "Not DEG", ], alpha = 0.7, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = c(
    "DEG in both" = "#2166AC",
    "DEG in softImpute only" = "#D6604D",
    "DEG in intersection only" = "#F4A582",
    "Not DEG" = "grey80")) +
  annotate("text", x = min(comp$logFC_intersect, na.rm = TRUE),
           y = max(comp$logFC_softimpute, na.rm = TRUE),
           label = sprintf("r = %.4f\nn = %s", r_val, formatC(nrow(comp), big.mark = ",")),
           hjust = 0, vjust = 1, size = 4) +
  labs(x = "logFC (intersection-only, 8,260 genes)",
       y = "logFC (softImpute, same 8,260 genes)",
       colour = "DEG status",
       title = "logFC agreement: softImpute vs intersection-only (shared genes)") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "fig_logfc_scatter_internal.pdf"), p_scatter, width = 8, height = 7)
png(file.path(output_dir, "fig_logfc_scatter_internal.png"),
    width = 2400, height = 2100, res = 300, type = "cairo")
print(p_scatter)
invisible(dev.off())
cat("Saved logFC scatter plot\n")

# --- Volcano plot: intersection-only full limma with softImpute-only DEGs highlighted ---

softimpute_only_degs <- setdiff(sig_softimpute$gene, sig_intersect$gene)
softimpute_only_in_shared <- softimpute_only_degs[softimpute_only_degs %in% shared_genes]

volcano_df <- data.frame(
  gene = full_intersect$gene,
  logFC = full_intersect$logFC,
  neg_log10_fdr = -log10(pmax(full_intersect$adj.P.Val, 1e-50)),
  stringsAsFactors = FALSE
)
volcano_df$neg_log10_fdr <- pmin(volcano_df$neg_log10_fdr, 14)
volcano_df$group <- ifelse(volcano_df$gene %in% softimpute_only_in_shared,
                           "DEG only in softImpute", "Other")
volcano_df <- volcano_df[order(volcano_df$group == "DEG only in softImpute"), ]

library(org.Hs.eg.db)
symbols <- mapIds(org.Hs.eg.db, keys = as.character(softimpute_only_in_shared),
                  column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
volcano_df$symbol <- symbols[volcano_df$gene]

top_highlight <- volcano_df[volcano_df$group == "DEG only in softImpute", ]
top_highlight <- top_highlight[order(-top_highlight$neg_log10_fdr), ][1:min(10, nrow(top_highlight)), ]

p_volcano <- ggplot(volcano_df, aes(x = logFC, y = neg_log10_fdr, colour = group)) +
  geom_point(data = volcano_df[volcano_df$group == "Other", ], alpha = 0.15, size = 0.8) +
  geom_point(data = volcano_df[volcano_df$group != "Other", ], alpha = 0.8, size = 2) +
  geom_text_repel(data = top_highlight, aes(label = symbol), size = 3,
                  max.overlaps = 15, show.legend = FALSE) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = c("Other" = "grey70", "DEG only in softImpute" = "#D6604D")) +
  labs(x = expression(log[2]~"fold change (intersection-only analysis)"),
       y = expression(-log[10]~"FDR"),
       colour = NULL,
       title = sprintf("Genes that become DEG only with softImpute (%d of %d shared genes)",
                        length(softimpute_only_in_shared), length(shared_genes))) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "fig_volcano_internal.pdf"), p_volcano, width = 8, height = 7)
png(file.path(output_dir, "fig_volcano_internal.png"),
    width = 2400, height = 2100, res = 300, type = "cairo")
print(p_volcano)
invisible(dev.off())
cat("Saved volcano plot\n")

# --- Clustered arrows: intersection-only → softImpute for DEGs in both ---

both_deg_ids <- intersect(sig_softimpute$gene, sig_intersect$gene)
both_deg_ids <- both_deg_ids[both_deg_ids %in% shared_genes]

fdr_cap <- 14

idx_io2 <- match(both_deg_ids, full_intersect$gene)
idx_si2 <- match(both_deg_ids, full_softimpute$gene)

arrows_df <- data.frame(
  gene = both_deg_ids,
  x0 = full_intersect$logFC[idx_io2],
  y0 = pmin(-log10(pmax(full_intersect$adj.P.Val[idx_io2], 1e-50)), fdr_cap),
  x1 = full_softimpute$logFC[idx_si2],
  y1 = pmin(-log10(pmax(full_softimpute$adj.P.Val[idx_si2], 1e-50)), fdr_cap),
  stringsAsFactors = FALSE
)

if (nrow(arrows_df) >= 8) {
  feat <- arrows_df[, c("x0", "y0", "x1", "y1")]
  feat_scaled <- scale(feat)

  best_k <- min(8, nrow(arrows_df) %/% 3)
  set.seed(42)
  km <- kmeans(feat_scaled, centers = best_k, nstart = 50, iter.max = 100)
  arrows_df$cluster <- km$cluster

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
  cat("\nClustered arrows: ", n_clusters, "clusters from", nrow(arrows_df), "shared DEGs\n")

  cluster_info <- do.call(rbind, lapply(sort(unique(arrows_df$cluster)), function(cl) {
    sub <- arrows_df[arrows_df$cluster == cl, ]
    pad_x <- 0.08; pad_y <- 0.15
    data.frame(
      cluster = cl, n = nrow(sub),
      from_cx = mean(sub$x0), from_cy = mean(sub$y0),
      to_cx = mean(sub$x1), to_cy = mean(sub$y1),
      from_xmin = min(sub$x0) - pad_x, from_xmax = max(sub$x0) + pad_x,
      from_ymin = min(sub$y0) - pad_y, from_ymax = max(sub$y0) + pad_y,
      to_xmin = min(sub$x1) - pad_x, to_xmax = max(sub$x1) + pad_x,
      to_ymin = min(sub$y1) - pad_y, to_ymax = max(sub$y1) + pad_y,
      stringsAsFactors = FALSE)
  }))

  bg_df <- data.frame(
    logFC = full_intersect$logFC,
    nlog10_fdr = pmin(-log10(pmax(full_intersect$adj.P.Val, 1e-50)), fdr_cap),
    stringsAsFactors = FALSE)

  pal <- c("#2166AC", "#D6604D", "#5AAE61", "#F4A582", "#8073AC",
           "#FDB863", "#1B9E77", "#E7298A", "#66A61E", "#7570B3")
  cluster_colours <- pal[seq_len(n_clusters)]
  names(cluster_colours) <- as.character(sort(unique(arrows_df$cluster)))

  p_arrows <- ggplot() +
    geom_point(data = bg_df, aes(x = logFC, y = nlog10_fdr),
               colour = "grey90", size = 0.2, alpha = 0.2) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey60", linewidth = 0.3) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey60", linewidth = 0.3)

  for (i in seq_len(nrow(cluster_info))) {
    ci <- cluster_info[i, ]
    col <- cluster_colours[as.character(ci$cluster)]
    p_arrows <- p_arrows +
      annotate("rect", xmin = ci$from_xmin, xmax = ci$from_xmax,
               ymin = ci$from_ymin, ymax = ci$from_ymax,
               fill = col, alpha = 0.10, colour = col, linetype = "dashed", linewidth = 0.5) +
      annotate("rect", xmin = ci$to_xmin, xmax = ci$to_xmax,
               ymin = ci$to_ymin, ymax = ci$to_ymax,
               fill = col, alpha = 0.15, colour = col, linetype = "solid", linewidth = 0.6) +
      annotate("segment", x = ci$from_cx, y = ci$from_cy,
               xend = ci$to_cx, yend = ci$to_cy,
               colour = col, linewidth = 1.2, alpha = 0.8,
               arrow = arrow(length = unit(3, "mm"), type = "closed")) +
      annotate("label", x = (ci$from_cx + ci$to_cx) / 2, y = (ci$from_cy + ci$to_cy) / 2,
               label = sprintf("n=%d", ci$n), size = 3, fontface = "bold", colour = col,
               fill = "white", alpha = 0.85, label.size = 0.3, label.padding = unit(1.5, "pt"))
  }

  p_arrows <- p_arrows +
    labs(x = expression(log[2]~"fold change"),
         y = expression(-log[10]~"FDR"),
         title = "Clustered DEG transitions: intersection-only → softImpute (shared genes)",
         subtitle = "Dashed = intersection-only, solid = softImpute") +
    coord_cartesian(xlim = range(c(bg_df$logFC, arrows_df$x0, arrows_df$x1)) * c(1.1, 1.1),
                    ylim = c(0, fdr_cap + 0.5)) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none", plot.title = element_text(size = 11))

  ggsave(file.path(output_dir, "fig_volcano_arrows_clustered_internal.pdf"), p_arrows, width = 11, height = 9)
  png(file.path(output_dir, "fig_volcano_arrows_clustered_internal.png"),
      width = 3300, height = 2700, res = 300, type = "cairo")
  print(p_arrows)
  invisible(dev.off())
  cat("Saved clustered arrows plot\n")
} else {
  cat("Too few shared DEGs for clustering (", nrow(arrows_df), ")\n")
}

# --- Summary ---

cat(sprintf("\n=== Internal concordance summary ===
Shared genes: %d
softImpute DEGs (total): %d
Intersection DEGs (total): %d
DEGs in both (shared genes): %d
DEGs only in softImpute (shared): %d
DEGs only in intersection (shared): %d
logFC Pearson r (shared): %.4f
Same direction (shared): %d/%d (%.1f%%)
\n",
  length(shared_genes),
  nrow(sig_softimpute), nrow(sig_intersect),
  length(both_deg_ids),
  length(softimpute_only_in_shared),
  length(setdiff(sig_intersect$gene, sig_softimpute$gene)),
  r_val,
  sum(comp$same_direction), nrow(comp),
  100 * sum(comp$same_direction) / nrow(comp)))

cat("All outputs saved to:", output_dir, "\n")
