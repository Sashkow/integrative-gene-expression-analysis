#!/usr/bin/env Rscript
#
# Recreate difexp_final.csv for 1_2 and 2_3 trimester comparisons.
# Run from the article_4_trims_compare_restoration/ directory.
#
# Usage:
#   cd scripts/integrative-gene-expression-analysis/scripts/article_4_trims_compare_restoration
#   Rscript run.R

library(igraph)
library(writexl)
library(dplyr)
library(RCurl)

# Source helper scripts from local scripts/ directory
base_dir <- getwd()
source(file.path(base_dir, "scripts", "centrality.R"))
source(file.path(base_dir, "scripts", "gene_metadata.R"))
source(file.path(base_dir, "scripts", "r_utils.R"))

# Load STRINGdb for its dependencies, then detach and replace with custom class
library(STRINGdb)
detach("package:STRINGdb", unload = TRUE)
source(file.path(base_dir, "scripts", "rstring.R"))


pipeline_restored <- function(analysis_id, confidence = 100) {
  cat("\n=== Running pipeline for", analysis_id, "===\n")

  data_path <- file.path(base_dir, "data", analysis_id)
  results_path <- file.path(base_dir, "results", analysis_id)

  if (!dir.exists(results_path)) {
    dir.create(results_path, recursive = TRUE)
  }

  setwd(results_path)

  # Read background expression
  sub_exprs <- read.csv(file.path(data_path, "sub_exprs.tsv"),
                        row.names = 1, sep = "\t")
  cat("Background genes:", nrow(sub_exprs), "\n")

  # Read DE results
  if (file.exists(file.path(data_path, "difexp.tsv"))) {
    difexp <- read.csv(file.path(data_path, "difexp.tsv"), sep = "\t")
  } else {
    difexp <- read.csv(file.path(data_path, "difexp.csv"))
    difexp <- difexp[!duplicated(difexp$SYMBOL), ]
  }
  rownames(difexp) <- difexp$SYMBOL
  cat("DE genes:", nrow(difexp), "\n")

  # Filter by logFC
  difexp <- difexp[which(abs(difexp$logFC) >= 1), ]
  cat("DE genes |logFC| >= 1:", nrow(difexp), "\n")

  # STRING database (cache downloads locally to avoid re-downloading)
  cache_dir <- file.path(base_dir, "stringdb_cache")
  if (!dir.exists(cache_dir)) dir.create(cache_dir)
  string_db <- STRINGdb$new(
    version = "11", species = 9606,
    score_threshold = confidence, input_directory = cache_dir
  )

  # Map background
  if ("X" %in% colnames(sub_exprs)) {
    sub_exprs$SYMBOL <- sub_exprs$X
  } else {
    sub_exprs$SYMBOL <- rownames(sub_exprs)
  }
  sub_exprs_mapped <- string_db$map(
    sub_exprs, "SYMBOL", removeUnmappedRows = TRUE
  )
  backgroundV <- sub_exprs_mapped$STRING_id
  string_db$set_background(backgroundV)
  cat("Background mapped to STRING:", nrow(sub_exprs_mapped), "\n")

  # Map DE genes (by SYMBOL)
  difexp_mapped <- string_db$map(
    difexp, "SYMBOL", removeUnmappedRows = TRUE, takeFirst = FALSE
  )
  difexp_mapped <- difexp_mapped[!duplicated(difexp_mapped$STRING_id), ]

  # Save unmapped genes
  unmapped_symbols <- setdiff(difexp$SYMBOL, difexp_mapped$SYMBOL)
  cat("Unmapped DE genes:", length(unmapped_symbols),
      "of", nrow(difexp), "\n")

  if (length(unmapped_symbols) > 0) {
    # Select available columns for unmapped report
    info_cols <- intersect(
      c("SYMBOL", "ENTREZID", "fENTREZID", "GENENAME"),
      colnames(difexp)
    )
    unmapped_info <- difexp[unmapped_symbols, info_cols, drop = FALSE]
    tryCatch({
      if (requireNamespace("biomaRt", quietly = TRUE)) {
        ensembl <- biomaRt::useEnsembl(
          "genes", dataset = "hsapiens_gene_ensembl"
        )
        bm <- biomaRt::getBM(
          attributes = c("hgnc_symbol", "gene_biotype"),
          filters = "hgnc_symbol",
          values = unmapped_info$SYMBOL, mart = ensembl
        )
        unmapped_info <- merge(
          unmapped_info, bm,
          by.x = "SYMBOL", by.y = "hgnc_symbol", all.x = TRUE
        )
        n_pc <- sum(
          unmapped_info$gene_biotype == "protein_coding",
          na.rm = TRUE
        )
        cat("  Protein-coding:", n_pc,
            " Non-protein-coding/unknown:",
            nrow(unmapped_info) - n_pc, "\n")
      }
    }, error = function(e) {
      cat("  (biomaRt unavailable:", conditionMessage(e), ")\n")
    })
    write.csv(unmapped_info,
              file.path(results_path, "unmapped_genes.csv"),
              row.names = FALSE)
    cat("  Saved unmapped_genes.csv\n")
  }
  cat("DE genes mapped to STRING:", nrow(difexp_mapped), "\n")

  # Build PPI subnetwork, remove isolated nodes
  G <- simplify(string_db$get_subnetwork(difexp_mapped$STRING_id))
  Isolated <- which(degree(G) == 0)
  G <- delete.vertices(G, Isolated)
  cat("Nodes on graph:", length(V(G)), "\n")

  # on_graph flag
  difexp_mapped$on_graph <- ifelse(
    difexp_mapped$STRING_id %in% names(V(G)), 1, 0
  )

  # Keep only on-graph genes
  difexp_mapped <- difexp_mapped[difexp_mapped$on_graph == 1, ]

  # Fastgreedy community detection -> cluster assignment
  sub_g <- G
  sub_fgreedy <- fastgreedy.community(
    sub_g, merges = TRUE, modularity = TRUE
  )
  cat("Modularity:", modularity(sub_fgreedy), "\n")

  cluster_membership <- membership(sub_fgreedy)
  difexp_mapped$cluster <- cluster_membership[difexp_mapped$STRING_id]

  # Up/down direction
  difexp_mapped$updown <- ifelse(
    difexp_mapped$logFC > 0, "green", "red"
  )

  # Betweenness centrality
  difexp_mapped <- add_cluster_max_btw_column(difexp_mapped, sub_g, string_db)

  # is_secreted
  difexp_mapped <- add_is_secreted_placenta_column(
    difexp_mapped, results_path
  )

  cat("Final genes:", nrow(difexp_mapped), "\n")

  # Write difexp_final.csv
  write.table(difexp_mapped, "difexp_final.csv",
              row.names = FALSE, sep = ",")

  cat("Done:", analysis_id, "\n\n")
  setwd(base_dir)
}


# Run both analyses
pipeline_restored("1_2")
pipeline_restored("2_3")

cat("\n=== All analyses complete ===\n")
cat("Results in results/1_2/ and results/2_3/\n")
cat("Compare with original_results/ to verify.\n")
