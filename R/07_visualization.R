#' Visualization Module
#'
#' Functions for network and cluster visualization
#'
#' @author Expression Integration Pipeline
#' @date 2025

#' Plot network with cluster colors
#'
#' @param G igraph network
#' @param difexp Differential expression with cluster and logFC
#' @param fgreedy Fastgreedy communities
#' @param output_file Output file path
#' @param layout_algo Layout algorithm (default: "fr")
#' @export
plot_network <- function(G, difexp, fgreedy, output_file,
                         layout_algo = "fr") {

  library(igraph)

  # Ensure genes are ordered to match network nodes
  difexp_ordered <- difexp[match(V(G)$name, difexp$STRING_id), ]

  # Set layout
  if (layout_algo == "fr") {
    G$layout <- layout_with_fr(G)
  } else if (layout_algo == "kk") {
    G$layout <- layout_with_kk(G)
  } else {
    G$layout <- layout_nicely(G)
  }

  # Node colors by cluster
  V(G)$color <- membership(fgreedy)

  # Node size by log fold-change
  if ("logFC" %in% colnames(difexp_ordered)) {
    V(G)$size <- abs(difexp_ordered$logFC) * 3 + 5
  } else {
    V(G)$size <- 8
  }

  # Edge weights
  E(G)$width <- E(G)$weight / 200

  # Plot
  png(output_file, width = 1200, height = 1200)
  plot(G,
       vertex.label = ifelse(V(G)$size > 15,
                             difexp_ordered$SYMBOL, ""),
       vertex.label.cex = 0.7,
       vertex.label.color = "black",
       edge.color = "gray80",
       main = "Protein-Protein Interaction Network")

  # Add legend
  legend("topright",
         legend = paste("Cluster", 1:length(fgreedy)),
         col = 1:length(fgreedy),
         pch = 19,
         cex = 0.8,
         title = "Communities")

  dev.off()

  message("Saved network plot: ", output_file)
  invisible(NULL)
}


#' Create volcano plot
#'
#' @param difexp Differential expression results
#' @param output_file Output file path
#' @param logfc_threshold Log fold-change threshold (default: 1)
#' @param pval_threshold P-value threshold (default: 0.05)
#' @param label_top Number of top genes to label (default: 20)
#' @export
plot_volcano <- function(difexp, output_file,
                         logfc_threshold = 1,
                         pval_threshold = 0.05,
                         label_top = 20) {

  library(ggplot2)

  if (!all(c("logFC", "adj.P.Val") %in% colnames(difexp))) {
    warning("Missing required columns for volcano plot")
    return(invisible(NULL))
  }

  # Add significance category
  difexp$significance <- "NS"
  difexp$significance[abs(difexp$logFC) > logfc_threshold] <- "LogFC"
  difexp$significance[difexp$adj.P.Val < pval_threshold] <- "P-value"
  difexp$significance[abs(difexp$logFC) > logfc_threshold &
                      difexp$adj.P.Val < pval_threshold] <- "Significant"

  # Get top genes to label
  difexp_sorted <- difexp[order(difexp$adj.P.Val), ]
  top_genes <- head(difexp_sorted, label_top)

  p <- ggplot(difexp, aes(x = logFC, y = -log10(adj.P.Val), color = significance)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = c("NS" = "gray", "LogFC" = "blue",
                                   "P-value" = "orange", "Significant" = "red")) +
    geom_vline(xintercept = c(-logfc_threshold, logfc_threshold),
               linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(pval_threshold),
               linetype = "dashed", color = "black") +
    theme_bw() +
    labs(title = "Volcano Plot",
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-value")

  # Add labels for top genes
  if ("SYMBOL" %in% colnames(top_genes)) {
    library(ggrepel)
    p <- p + geom_text_repel(data = top_genes,
                             aes(label = SYMBOL),
                             size = 3, max.overlaps = 20)
  }

  ggsave(output_file, p, width = 10, height = 8)
  message("Saved volcano plot: ", output_file)

  invisible(p)
}


#' Create cluster size barplot
#'
#' @param summary_data Cluster summary data frame
#' @param output_file Output file path
#' @export
plot_cluster_sizes <- function(summary_data, output_file) {

  library(ggplot2)

  # Check if we have data
  if (is.null(summary_data) || nrow(summary_data) == 0) {
    message("No cluster data to plot")
    return(invisible(NULL))
  }

  p <- ggplot(summary_data, aes(x = factor(cluster), y = genes_total)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_text(aes(label = genes_total), vjust = -0.5) +
    theme_bw() +
    labs(title = "Cluster Sizes",
         x = "Cluster",
         y = "Number of Genes")

  # Save with explicit device to avoid graphics API issues
  tryCatch({
    ggsave(output_file, p, width = 10, height = 6, device = "png")
    message("Saved cluster size plot: ", output_file)
  }, error = function(e) {
    warning("Failed to save cluster size plot: ", e$message)
  })

  invisible(p)
}


#' Create enrichment coverage plot
#'
#' @param summary_data Cluster summary data frame
#' @param output_file Output file path
#' @export
plot_enrichment_coverage <- function(summary_data, output_file) {

  library(ggplot2)

  # Check if we have data
  if (is.null(summary_data) || nrow(summary_data) == 0) {
    message("No enrichment data to plot")
    return(invisible(NULL))
  }

  p <- ggplot(summary_data, aes(x = factor(cluster), y = coverage)) +
    geom_bar(stat = "identity", fill = "coral") +
    geom_text(aes(label = paste0(round(coverage * 100, 1), "%")),
              vjust = -0.5) +
    theme_bw() +
    labs(title = "Enrichment Coverage by Cluster",
         x = "Cluster",
         y = "Proportion of Genes Covered") +
    ylim(0, 1)

  # Save with explicit device to avoid graphics API issues
  tryCatch({
    ggsave(output_file, p, width = 10, height = 6, device = "png")
    message("Saved enrichment coverage plot: ", output_file)
  }, error = function(e) {
    warning("Failed to save enrichment coverage plot: ", e$message)
  })

  invisible(p)
}


#' Create heatmap of top differential genes
#'
#' @param exprs Expression matrix
#' @param difexp Differential expression results
#' @param pdata Phenodata
#' @param group_col Column for grouping (default: "trim_term")
#' @param top_n Number of top genes (default: 50)
#' @param output_file Output file path
#' @export
plot_difexp_heatmap <- function(exprs, difexp, pdata, group_col = "trim_term",
                                 top_n = 50, output_file) {

  library(pheatmap)

  # Get top genes
  if ("adj.P.Val" %in% colnames(difexp)) {
    top_genes <- head(difexp[order(difexp$adj.P.Val), ], top_n)
  } else {
    top_genes <- head(difexp, top_n)
  }

  # Get ENTREZID for indexing (expression matrix is indexed by ENTREZID)
  if (!"ENTREZID" %in% colnames(top_genes)) {
    warning("No ENTREZID column found in differential expression results")
    return(invisible(NULL))
  }

  # Debug information
  message("Top genes ENTREZID (first 5): ", paste(head(top_genes$ENTREZID, 5), collapse=", "))
  message("Expression matrix rownames (first 5): ", paste(head(rownames(exprs), 5), collapse=", "))
  message("Expression matrix class: ", class(exprs))

  # Find genes that exist in expression matrix
  gene_ids <- intersect(as.character(top_genes$ENTREZID), rownames(exprs))

  if (length(gene_ids) == 0) {
    warning("No genes found in expression matrix")
    warning("  Top genes ENTREZIDs: ", paste(head(top_genes$ENTREZID, 10), collapse=", "))
    warning("  Expression rownames: ", paste(head(rownames(exprs), 10), collapse=", "))
    return(invisible(NULL))
  }

  message("Plotting heatmap for ", length(gene_ids), " out of ", nrow(top_genes), " top genes")

  # Subset expression matrix
  exprs_subset <- exprs[gene_ids, , drop = FALSE]

  # Use gene symbols as row names if available
  if ("SYMBOL" %in% colnames(top_genes)) {
    idx <- match(gene_ids, top_genes$ENTREZID)
    symbols <- top_genes$SYMBOL[idx]
    # Use symbol if available, otherwise ENTREZID
    rownames(exprs_subset) <- ifelse(!is.na(symbols), symbols, gene_ids)
  }

  # Create annotation - align with expression columns
  sample_names <- colnames(exprs_subset)
  pdata_subset <- pdata[make.names(pdata$arraydatafile_exprscolumnnames) %in% sample_names, ]

  annotation_col <- data.frame(
    Group = pdata_subset[[group_col]],
    row.names = make.names(pdata_subset$arraydatafile_exprscolumnnames)
  )

  # Plot
  png(output_file, width = 1200, height = 1000)
  pheatmap(exprs_subset,
           scale = "row",
           annotation_col = annotation_col,
           show_colnames = FALSE,
           main = paste("Top", top_n, "Differential Genes"))
  dev.off()

  message("Saved heatmap: ", output_file)
  invisible(NULL)
}
