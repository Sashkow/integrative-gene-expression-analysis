#' Network Clustering Module
#'
#' Functions for STRING database mapping and fastgreedy clustering
#'
#' @author Expression Integration Pipeline
#' @date 2025

#' Initialize STRING database
#'
#' @param version STRING version (default: "11")
#' @param species NCBI taxonomy ID (default: 9606 for human)
#' @param score_threshold Confidence score threshold (default: 400)
#' @param use_cache Use cached database if available (default: TRUE)
#' @param cache_dir Directory for cache (default: "data/stringdb_cache")
#' @return STRINGdb object
#' @export
initialize_stringdb <- function(version = "11", species = 9606, score_threshold = 400,
                                 use_cache = TRUE, cache_dir = "data/stringdb_cache") {

  library(STRINGdb)

  # Set up cache directory
  input_dir <- ""
  if (use_cache) {
    # Create cache directory if it doesn't exist
    if (!dir.exists(cache_dir)) {
      dir.create(cache_dir, recursive = TRUE)
      message("Created cache directory: ", cache_dir)
    }

    # Check if cache exists
    cache_file <- file.path(cache_dir,
                            paste0(species, ".protein.links.v", version, ".txt.gz"))
    if (file.exists(cache_file)) {
      message("Using cached STRING database from: ", cache_dir)
    } else {
      message("Cache not found - will download to: ", cache_dir)
    }
    input_dir <- cache_dir
  } else {
    message("Cache disabled - using temporary directory")
  }

  string_db <- STRINGdb$new(
    version = version,
    species = species,
    score_threshold = score_threshold,
    input_directory = input_dir
  )

  message("Initialized STRING database:")
  message("  Version: ", version)
  message("  Species: ", species)
  message("  Score threshold: ", score_threshold)

  return(string_db)
}


#' Map genes to STRING database
#'
#' @param difexp Differential expression results
#' @param background_exprs Background expression matrix for enrichment
#' @param string_db STRINGdb object
#' @param gene_col Column containing gene symbols (default: "SYMBOL")
#' @return List with mapped difexp and background
#' @export
map_to_stringdb <- function(difexp, background_exprs, string_db,
                             gene_col = "SYMBOL") {

  # Map differential genes
  difexp_mapped <- string_db$map(difexp, gene_col,
                                  removeUnmappedRows = TRUE,
                                  takeFirst = FALSE)
  difexp_mapped <- difexp_mapped[!duplicated(difexp_mapped$STRING_id), ]

  message("Mapped differential genes: ", nrow(difexp), " -> ",
          nrow(difexp_mapped), " STRING IDs")

  # Map background
  if (!is.null(background_exprs)) {
    # Convert ENTREZID to SYMBOL for background
    library(org.Hs.eg.db)
    library(AnnotationDbi)

    entrez_ids <- rownames(background_exprs)
    symbols <- AnnotationDbi::select(org.Hs.eg.db,
                                      keys = entrez_ids,
                                      columns = "SYMBOL",
                                      keytype = "ENTREZID")
    symbols <- symbols[!duplicated(symbols$ENTREZID) & !is.na(symbols$SYMBOL), ]

    background_df <- data.frame(SYMBOL = symbols$SYMBOL)
    background_mapped <- string_db$map(background_df, "SYMBOL",
                                        removeUnmappedRows = TRUE)
    background_ids <- unique(background_mapped$STRING_id)

    string_db$set_background(background_ids)
    message("Set background: ", length(background_ids), " genes")
  }

  return(list(
    difexp = difexp_mapped,
    background = if(!is.null(background_exprs)) background_ids else NULL
  ))
}


#' Build protein-protein interaction network
#'
#' @param string_ids Vector of STRING IDs
#' @param string_db STRINGdb object
#' @return igraph object
#' @export
build_ppi_network <- function(string_ids, string_db) {

  library(igraph)

  # Check if we have any genes
  if (length(string_ids) == 0) {
    warning("No STRING IDs provided - returning empty network")
    return(make_empty_graph(n = 0))
  }

  # Get subnetwork
  G <- string_db$get_subnetwork(string_ids)
  G <- simplify(G)

  message("Built PPI network: ", vcount(G), " nodes, ", ecount(G), " edges")

  if (vcount(G) == 0) {
    warning("No interactions found in STRING database - network is empty")
    warning("  Try lowering score_threshold in config")
    return(G)
  }

  # Remove isolated nodes
  isolated <- which(degree(G) == 0)
  if (length(isolated) > 0) {
    G <- delete.vertices(G, isolated)
    message("  Removed ", length(isolated), " isolated nodes")
    message("  Final network: ", vcount(G), " nodes, ", ecount(G), " edges")
  }

  return(G)
}


#' Perform fastgreedy community detection
#'
#' @param G igraph network object
#' @return fastgreedy communities object
#' @export
fastgreedy_clustering <- function(G) {

  library(igraph)

  # Check if network is empty
  if (vcount(G) == 0) {
    warning("Empty network - cannot perform clustering")
    return(NULL)
  }

  fgreedy <- fastgreedy.community(G, merges = TRUE, modularity = TRUE)

  message("Fastgreedy clustering:")
  message("  Communities: ", length(fgreedy))
  message("  Modularity: ", round(modularity(fgreedy), 3))
  message("  Sizes: ", paste(sizes(fgreedy), collapse = ", "))

  return(fgreedy)
}


#' Add cluster assignments to differential expression results
#'
#' @param difexp Mapped differential expression results with STRING_id
#' @param fgreedy Fastgreedy communities object
#' @param G igraph network
#' @return difexp with cluster column
#' @export
add_cluster_column <- function(difexp, fgreedy, G) {

  library(igraph)

  # Check if network is empty or clustering failed
  if (is.null(fgreedy) || vcount(G) == 0) {
    warning("Empty network or no clustering - cannot assign clusters")
    difexp$cluster <- NA
    difexp$on_graph <- 0
    return(difexp)
  }

  # Get membership
  com <- data.frame(
    STRING_id = V(G)$name,
    cluster = membership(fgreedy),
    stringsAsFactors = FALSE
  )

  # Merge with difexp
  difexp_clustered <- merge(difexp, com, by = "STRING_id", all.x = TRUE)
  difexp_clustered$cluster <- as.numeric(difexp_clustered$cluster)

  # Mark genes on graph
  difexp_clustered$on_graph <- ifelse(!is.na(difexp_clustered$cluster), 1, 0)

  message("Cluster assignment:")
  message("  Genes on graph: ", sum(difexp_clustered$on_graph))
  message("  Clusters: ", table(difexp_clustered$cluster))

  return(difexp_clustered)
}


#' Recursively cluster each community
#'
#' @param difexp Differential expression with clusters
#' @param string_db STRINGdb object
#' @param output_dir Output directory
#' @param min_cluster_size Minimum cluster size (default: 3)
#' @export
cluster_recursive <- function(difexp, string_db, output_dir,
                               min_cluster_size = 3) {

  library(igraph)

  clusters_dir <- file.path(output_dir, "clusters")
  if (!dir.exists(clusters_dir)) {
    dir.create(clusters_dir, recursive = TRUE)
  }

  # Get unique clusters
  cluster_ids <- sort(unique(difexp$cluster[!is.na(difexp$cluster)]))

  for (cluster_id in cluster_ids) {

    cluster_genes <- difexp[difexp$cluster == cluster_id, ]

    if (nrow(cluster_genes) < min_cluster_size) {
      message("Skipping cluster ", cluster_id, " (too small: ", nrow(cluster_genes), ")")
      next
    }

    message("\nProcessing cluster ", cluster_id, " (", nrow(cluster_genes), " genes)")

    # Create cluster directory
    cluster_dir <- file.path(clusters_dir, cluster_id)
    dir.create(cluster_dir, showWarnings = FALSE)

    # Build subnetwork
    sub_g <- build_ppi_network(cluster_genes$STRING_id, string_db)

    if (vcount(sub_g) >= min_cluster_size) {
      # Cluster subnetwork
      sub_fgreedy <- fastgreedy_clustering(sub_g)

      # Add sub-cluster assignments (rename cluster column to avoid conflict)
      if (!is.null(sub_fgreedy)) {
        # Remove old cluster column temporarily
        parent_cluster <- cluster_genes$cluster
        cluster_genes$cluster <- NULL

        # Add new sub-cluster assignments
        cluster_genes <- add_cluster_column(cluster_genes, sub_fgreedy, sub_g)

        # Rename to subcluster
        colnames(cluster_genes)[colnames(cluster_genes) == "cluster"] <- "subcluster"
        cluster_genes$parent_cluster <- parent_cluster
      }

      # Save cluster data
      write.csv(cluster_genes,
                file.path(cluster_dir, "cluster_genes.csv"),
                row.names = FALSE)
    } else {
      message("  Subnetwork too small for clustering")
      # Add subcluster column as NA
      cluster_genes$subcluster <- NA
      cluster_genes$parent_cluster <- cluster_genes$cluster

      write.csv(cluster_genes,
                file.path(cluster_dir, "cluster_genes.csv"),
                row.names = FALSE)
    }
  }

  message("\nRecursive clustering completed")
  invisible(NULL)
}
