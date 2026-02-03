#' Enrichment Analysis Module
#'
#' Functions for functional enrichment analysis
#'
#' @author Expression Integration Pipeline
#' @date 2025

#' Perform STRING enrichment analysis
#'
#' @param string_ids Vector of STRING IDs
#' @param string_db STRINGdb object
#' @param max_pvalue P-value threshold (default: 0.05)
#' @param category Enrichment category (default: "Process")
#' @return Data frame with enrichment results
#' @export
string_enrichment <- function(string_ids, string_db,
                               max_pvalue = 0.05,
                               category = "Process") {

  # Check if we have any genes
  if (length(string_ids) == 0) {
    message("No STRING IDs provided for enrichment")
    return(data.frame())
  }

  # Try enrichment analysis
  enrichment_all <- tryCatch({
    string_db$get_enrichment(string_ids)
  }, error = function(e) {
    warning("Enrichment analysis failed: ", e$message)
    return(data.frame())
  })

  if (is.null(enrichment_all) || nrow(enrichment_all) == 0) {
    message("No enrichment results found")
    return(data.frame())
  }

  # Filter by category and p-value
  enrichment <- enrichment_all[enrichment_all$category == category, ]
  enrichment <- enrichment[enrichment$p_value < max_pvalue, ]

  # Format for REVIGO
  enrichment$revigo_term <- gsub("\\.", ":", enrichment$term)

  # Parse gene lists
  if (nrow(enrichment) > 0) {
    if (grepl(",", enrichment$inputGenes[1])) {
      enrichment$list_genes <- strsplit(enrichment$inputGenes, ",")
    } else {
      enrichment$list_genes <- strsplit(enrichment$inputGenes, ";")
    }
  }

  message("Enrichment results: ", nrow(enrichment), " terms (", category,
          ", p < ", max_pvalue, ")")

  return(enrichment)
}


#' Calculate enrichment coverage
#'
#' Calculate what proportion of genes are covered by enrichment terms
#'
#' @param difexp Differential expression results with STRING_id
#' @param enrichment Enrichment results with list_genes column
#' @return List with coverage metrics
#' @export
calculate_enrichment_coverage <- function(difexp, enrichment) {

  if (nrow(enrichment) == 0) {
    return(list(
      coverage = 0,
      genes_covered = 0,
      genes_total = nrow(difexp)
    ))
  }

  # Get union of all genes in enrichment
  gene_union <- unique(unlist(enrichment$list_genes))

  # Calculate coverage
  genes_in_difexp <- sum(difexp$STRING_id %in% gene_union)
  coverage <- genes_in_difexp / nrow(difexp)

  message("Enrichment coverage: ", genes_in_difexp, "/", nrow(difexp),
          " genes (", round(coverage * 100, 1), "%)")

  return(list(
    coverage = coverage,
    genes_covered = genes_in_difexp,
    genes_total = nrow(difexp),
    gene_union = gene_union
  ))
}


#' Cluster-by-cluster enrichment analysis
#'
#' @param difexp Differential expression with cluster assignments
#' @param string_db STRINGdb object
#' @param output_dir Output directory
#' @param max_pvalue P-value threshold (default: 0.05)
#' @export
cluster_enrichment <- function(difexp, string_db, output_dir,
                                max_pvalue = 0.05) {

  clusters_dir <- file.path(output_dir, "clusters")
  summary_data <- data.frame()

  cluster_ids <- sort(unique(difexp$cluster[!is.na(difexp$cluster)]))

  for (cluster_id in cluster_ids) {

    cluster_genes <- difexp[difexp$cluster == cluster_id & !is.na(difexp$cluster), ]

    message("\nEnriching cluster ", cluster_id, " (", nrow(cluster_genes), " genes)")

    cluster_dir <- file.path(clusters_dir, cluster_id)
    if (!dir.exists(cluster_dir)) {
      dir.create(cluster_dir, recursive = TRUE)
    }

    # Perform enrichment
    enrichment <- string_enrichment(
      cluster_genes$STRING_id,
      string_db,
      max_pvalue = max_pvalue
    )

    # Calculate coverage
    coverage_metrics <- calculate_enrichment_coverage(cluster_genes, enrichment)

    # Save enrichment results
    if (nrow(enrichment) > 0) {
      enrichment_save <- enrichment[, !names(enrichment) %in% c("list_genes")]
      write.csv(enrichment_save,
                file.path(cluster_dir, "cluster_enrichment.csv"),
                row.names = FALSE)

      # Save for REVIGO
      revigo_df <- enrichment[, c("revigo_term", "p_value")]
      write.table(revigo_df,
                  file.path(cluster_dir, "cluster_to_revigo.tsv"),
                  sep = " ", row.names = FALSE, quote = FALSE)
    }

    # Update summary
    summary_row <- data.frame(
      cluster = cluster_id,
      genes_total = nrow(cluster_genes),
      enrichment_terms = nrow(enrichment),
      coverage = coverage_metrics$coverage,
      genes_covered = coverage_metrics$genes_covered
    )

    summary_data <- rbind(summary_data, summary_row)
  }

  # Save summary
  write.csv(summary_data,
            file.path(output_dir, "enrichment_summary.csv"),
            row.names = FALSE)

  message("\nCluster enrichment completed")
  return(summary_data)
}


#' Add gene-term membership counts
#'
#' Calculate how many enrichment terms each gene belongs to
#'
#' @param difexp Differential expression results
#' @param enrichment Enrichment results
#' @return difexp with term_count column
#' @export
add_term_membership <- function(difexp, enrichment) {

  difexp$term_count <- 0

  if (is.null(enrichment) || nrow(enrichment) == 0) {
    message("No enrichment data - term_count set to 0")
    return(difexp)
  }

  if (!"list_genes" %in% colnames(enrichment)) {
    warning("No list_genes column in enrichment - term_count set to 0")
    return(difexp)
  }

  if (!"STRING_id" %in% colnames(difexp)) {
    warning("No STRING_id column in difexp - term_count set to 0")
    return(difexp)
  }

  for (i in 1:nrow(difexp)) {
    gene <- difexp$STRING_id[i]
    if (!is.na(gene)) {
      count <- sum(sapply(enrichment$list_genes, function(x) gene %in% x))
      difexp$term_count[i] <- count
    }
  }

  message("Added term membership counts")
  message("  Range: ", min(difexp$term_count), " - ", max(difexp$term_count))

  return(difexp)
}


#' Perform enrichR analysis
#'
#' @param gene_symbols Vector of gene symbols
#' @param databases Vector of enrichR database names
#' @param max_pvalue P-value threshold (default: 0.05)
#' @return List of enrichment results
#' @export
enrichr_analysis <- function(gene_symbols, databases, max_pvalue = 0.05) {

  library(enrichR)

  message("Running enrichR with databases: ", paste(databases, collapse = ", "))

  enriched_full <- enrichr(gene_symbols, databases)

  results <- list()
  for (db in databases) {
    if (db %in% names(enriched_full)) {
      enriched <- enriched_full[[db]]
      enriched <- enriched[enriched$Adjusted.P.value <= max_pvalue, ]

      message("  ", db, ": ", nrow(enriched), " terms")
      results[[db]] <- enriched
    }
  }

  return(results)
}
