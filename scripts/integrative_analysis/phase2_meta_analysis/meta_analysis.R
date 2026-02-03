#' Phase 2: Meta-Analysis Methods
#'
#' Wrapper functions for DExMA and RankProd meta-analysis methods.
#' Provides consistent interface for running meta-analysis on multiple studies.
#'
#' @author Expression Integration Pipeline
#' @date 2025

# ============================================================================
# Dependencies
# ============================================================================

suppressPackageStartupMessages({
  library(limma)
})

# ============================================================================
# Data Preparation Functions
# ============================================================================

#' Load datasets and phenodata for meta-analysis
#'
#' @param mapped_path Directory containing expression TSV files
#' @param phenodata_path Path to sample metadata CSV
#' @param datasets Vector of dataset IDs to include (or NULL for all)
#' @param group_column Column name for group variable
#' @param baseline Baseline group name
#' @param contrast Contrast group name
#' @return List with expression matrices and phenotype info per study
#' @export
load_studies_for_meta <- function(mapped_path,
                                   phenodata_path,
                                   datasets = NULL,
                                   group_column = "trimester",
                                   baseline = "Second_Trimester",
                                   contrast = "Third_Trimester") {

  # Load phenodata
  phenodata <- read.csv(phenodata_path, stringsAsFactors = FALSE)

  # Get expression files
  expr_files <- list.files(mapped_path, pattern = "\\.tsv$", full.names = TRUE)

  if (!is.null(datasets)) {
    expr_files <- expr_files[sapply(datasets, function(d) {
      any(grepl(d, expr_files))
    })]
  }

  studies <- list()

  for (f in expr_files) {
    # Extract dataset ID from filename
    dataset_id <- gsub(".*/(GSE[0-9]+).*\\.tsv$", "\\1", f)

    # Load expression data
    expr <- read.table(f, header = TRUE, sep = "\t", row.names = 1,
                       check.names = FALSE, stringsAsFactors = FALSE)

    # Match samples to phenodata
    sample_names <- colnames(expr)
    pd_matched <- phenodata[phenodata$arraydatafile_exprscolumnnames %in% sample_names, ]

    if (nrow(pd_matched) == 0) {
      warning("No matching phenodata for ", dataset_id, ", skipping")
      next
    }

    # Filter to samples with group info
    pd_matched <- pd_matched[!is.na(pd_matched[[group_column]]), ]
    pd_matched <- pd_matched[pd_matched[[group_column]] %in% c(baseline, contrast), ]

    if (nrow(pd_matched) < 3) {
      warning("Fewer than 3 samples for ", dataset_id, ", skipping")
      next
    }

    # Subset expression to matched samples
    matched_samples <- pd_matched$arraydatafile_exprscolumnnames
    expr <- expr[, colnames(expr) %in% matched_samples, drop = FALSE]

    # Reorder phenodata to match expression columns
    pd_matched <- pd_matched[match(colnames(expr), pd_matched$arraydatafile_exprscolumnnames), ]

    # Create group vector (0 = baseline, 1 = contrast)
    groups <- ifelse(pd_matched[[group_column]] == contrast, 1, 0)

    studies[[dataset_id]] <- list(
      expr = as.matrix(expr),
      pheno = pd_matched,
      groups = groups,
      n_baseline = sum(groups == 0),
      n_contrast = sum(groups == 1)
    )

    cat("  Loaded", dataset_id, ":",
        ncol(expr), "samples,",
        nrow(expr), "genes,",
        sum(groups == 0), baseline, "/",
        sum(groups == 1), contrast, "\n")
  }

  return(studies)
}


#' Filter studies to only those with samples in both groups
#'
#' Required for DExMA and RankProd which need balanced studies.
#'
#' @param studies List of study objects from load_studies_for_meta
#' @param min_per_group Minimum samples required per group (default: 2)
#' @return Filtered studies list
#' @export
filter_balanced_studies <- function(studies, min_per_group = 2) {
  balanced <- list()

  for (study_id in names(studies)) {
    study <- studies[[study_id]]
    if (study$n_baseline >= min_per_group && study$n_contrast >= min_per_group) {
      balanced[[study_id]] <- study
    }
  }

  n_removed <- length(studies) - length(balanced)
  if (n_removed > 0) {
    cat("Filtered to balanced studies:", length(balanced), "of", length(studies),
        "(removed", n_removed, "with <", min_per_group, "samples in one group)\n")
  }

  return(balanced)
}


#' Filter studies to protein-coding genes
#'
#' @param studies List of study objects from load_studies_for_meta
#' @return Filtered studies list
#' @export
filter_studies_protein_coding <- function(studies) {
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("org.Hs.eg.db package required for gene filtering")
  }

  # Get protein-coding gene IDs
  gene_types <- AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, "ENTREZID"),
    columns = c("ENTREZID", "GENETYPE"),
    keytype = "ENTREZID"
  )
  protein_coding <- gene_types$ENTREZID[gene_types$GENETYPE == "protein-coding"]

  total_before <- 0
  total_after <- 0

 for (study_id in names(studies)) {
    genes_before <- rownames(studies[[study_id]]$expr)
    total_before <- total_before + length(genes_before)

    genes_keep <- genes_before[genes_before %in% protein_coding]
    studies[[study_id]]$expr <- studies[[study_id]]$expr[genes_keep, , drop = FALSE]

    total_after <- total_after + length(genes_keep)
  }

  cat("Filtered to protein-coding genes:",
      total_before, "->", total_after, "total gene entries\n")

  return(studies)
}


# ============================================================================
# Per-Study Differential Expression
# ============================================================================

#' Run limma DE analysis on a single study
#'
#' @param expr Expression matrix (genes x samples)
#' @param groups Binary vector (0 = baseline, 1 = contrast)
#' @return Data frame with limma results
#' @export
run_limma_single_study <- function(expr, groups) {
  design <- model.matrix(~ groups)

  fit <- lmFit(expr, design)
  fit <- eBayes(fit)

  results <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  results$gene_id <- rownames(results)

  return(results)
}


#' Run limma on all studies
#'
#' @param studies List of study objects
#' @return List of DE results per study
#' @export
run_per_study_de <- function(studies) {
  results <- list()

  for (study_id in names(studies)) {
    cat("  Running limma on", study_id, "...\n")

    study <- studies[[study_id]]
    de <- run_limma_single_study(study$expr, study$groups)
    de$study <- study_id

    results[[study_id]] <- de
  }

  return(results)
}


# ============================================================================
# DExMA Meta-Analysis
# ============================================================================

#' Prepare data for DExMA
#'
#' @param studies List of study objects
#' @return Object suitable for DExMA createObjectMA
#' @export
prepare_dexma_input <- function(studies) {
  if (!requireNamespace("DExMA", quietly = TRUE)) {
    stop("DExMA package not installed. Run: BiocManager::install('DExMA')")
  }

  # DExMA expects a list of lists, each with:
  # - GEX: expression matrix
  # - Sample_Pheno: data frame with phenotype
  # - Var_Phenotype: name of the phenotype column

  dexma_list <- list()

  for (study_id in names(studies)) {
    study <- studies[[study_id]]

    # Create phenotype data frame
    pheno_df <- data.frame(
      sample = colnames(study$expr),
      group = factor(study$groups, levels = c(0, 1),
                     labels = c("control", "case"))
    )
    rownames(pheno_df) <- pheno_df$sample

    dexma_list[[study_id]] <- list(
      GEX = study$expr,
      Sample_Pheno = pheno_df,
      Var_Phenotype = "group"
    )
  }

  return(dexma_list)
}


#' Run DExMA meta-analysis
#'
#' @param studies List of study objects from load_studies_for_meta
#' @param effect_size "SMD" or "MD"
#' @param missAllow Maximum proportion of missing studies per gene
#' @param impute Whether to impute missing values
#' @return DExMA results object
#' @export
run_dexma_meta <- function(studies,
                           effect_size = "SMD",
                           missAllow = 0.3,
                           impute = TRUE) {

  if (!requireNamespace("DExMA", quietly = TRUE)) {
    stop("DExMA package not installed")
  }

  cat("\n=== Running DExMA Meta-Analysis ===\n\n")
  cat("Effect size:", effect_size, "\n")
  cat("Missing studies allowed:", missAllow * 100, "%\n")
  cat("Imputation:", ifelse(impute, "enabled", "disabled"), "\n\n")

  # Prepare input
  dexma_input <- prepare_dexma_input(studies)

  # Create meta-analysis object
  cat("Creating meta-analysis object...\n")

  # Extract expression matrices and phenotype info
  listEX <- lapply(dexma_input, function(x) x$GEX)
  listPheno <- lapply(dexma_input, function(x) x$Sample_Pheno)

  # Debug: print structure
  cat("  Studies:", length(listEX), "\n")
  for (i in seq_along(listEX)) {
    nm <- names(listEX)[i]
    cat("  ", nm, ": ", nrow(listEX[[i]]), " genes, ", ncol(listEX[[i]]),
        " samples, groups: ", paste(unique(listPheno[[i]]$group), collapse = "/"),
        "\n", sep = "")
  }

  # Create object - expGroups/refGroups are the actual group labels
  meta_obj <- DExMA::createObjectMA(
    listEX = listEX,
    listPheno = listPheno,
    namePheno = rep("group", length(listEX)),
    expGroups = rep("case", length(listEX)),
    refGroups = rep("control", length(listEX))
  )

  # Run meta-analysis
  cat("Running effect-size meta-analysis...\n")
  meta_results <- DExMA::metaAnalysisDE(
    meta_obj,
    typeMethod = "REM",       # Random effects model
    missAllow = missAllow,
    proportionData = 1 - missAllow
  )

  # Extract results (DExMA column names: Com.ES, ES.var, Pval, FDR, propDataset)
  results_df <- data.frame(
    gene_id = rownames(meta_results),
    logFC = meta_results$Com.ES,
    se = sqrt(meta_results$ES.var),
    pvalue = meta_results$Pval,
    fdr = meta_results$FDR,
    prop_datasets = meta_results$propDataset,
    stringsAsFactors = FALSE
  )

  results_df <- results_df[order(results_df$fdr), ]

  cat("\nDExMA complete.\n")
  cat("Genes analyzed:", nrow(results_df), "\n")
  cat("Significant (FDR < 0.05):", sum(results_df$fdr < 0.05, na.rm = TRUE), "\n")

  return(results_df)
}


# ============================================================================
# RankProd Meta-Analysis
# ============================================================================

#' Prepare data for RankProd
#'
#' @param studies List of study objects
#' @param min_studies Minimum studies a gene must appear in
#' @return List with combined matrix and class labels
#' @export
prepare_rankprod_input <- function(studies, min_studies = 4) {
  if (!requireNamespace("RankProd", quietly = TRUE)) {
    stop("RankProd package not installed. Run: BiocManager::install('RankProd')")
  }

  # Get all genes
  all_genes <- unique(unlist(lapply(studies, function(s) rownames(s$expr))))

  # Count studies per gene
  gene_counts <- table(unlist(lapply(studies, function(s) rownames(s$expr))))
  genes_keep <- names(gene_counts[gene_counts >= min_studies])

  cat("Genes in at least", min_studies, "studies:", length(genes_keep), "\n")

  # Build combined matrix (genes x samples) and origin vector
  combined_expr <- NULL
  combined_groups <- c()
  combined_origin <- c()

  for (i in seq_along(studies)) {
    study_id <- names(studies)[i]
    study <- studies[[study_id]]

    # Subset to common genes
    expr_sub <- study$expr[rownames(study$expr) %in% genes_keep, , drop = FALSE]

    if (is.null(combined_expr)) {
      combined_expr <- expr_sub
    } else {
      # Align genes
      common_genes <- intersect(rownames(combined_expr), rownames(expr_sub))
      combined_expr <- cbind(
        combined_expr[common_genes, , drop = FALSE],
        expr_sub[common_genes, , drop = FALSE]
      )
    }

    combined_groups <- c(combined_groups, study$groups)
    combined_origin <- c(combined_origin, rep(i, ncol(study$expr)))
  }

  return(list(
    expr = combined_expr,
    groups = combined_groups,
    origin = combined_origin,
    n_genes = nrow(combined_expr),
    n_samples = ncol(combined_expr)
  ))
}


#' Run RankProd meta-analysis
#'
#' @param studies List of study objects
#' @param min_studies Minimum studies a gene must appear in
#' @param num_perm Number of permutations
#' @param logged Is data already log-transformed?
#' @return RankProd results data frame
#' @export
run_rankprod_meta <- function(studies,
                              min_studies = 4,
                              num_perm = 1000,
                              logged = TRUE) {

  if (!requireNamespace("RankProd", quietly = TRUE)) {
    stop("RankProd package not installed")
  }

  cat("\n=== Running RankProd Meta-Analysis ===\n\n")
  cat("Minimum studies per gene:", min_studies, "\n")
  cat("Permutations:", num_perm, "\n")
  cat("Data logged:", logged, "\n\n")

  # Prepare input
  rp_input <- prepare_rankprod_input(studies, min_studies)

  cat("Combined matrix:", rp_input$n_genes, "genes x", rp_input$n_samples, "samples\n")
  cat("Studies:", length(unique(rp_input$origin)), "\n\n")

  # Run RankProd (use RPadvance for multi-study analysis)
  cat("Running RPadvance (this may take a while)...\n")
  rp_result <- RankProd::RPadvance(
    data = rp_input$expr,
    cl = rp_input$groups,
    origin = rp_input$origin,
    num.perm = num_perm,
    logged = logged,
    na.rm = TRUE,
    plot = FALSE
  )

  # Extract results
  # RankProd returns separate results for up and down regulation
  # pfp = percentage of false positives (similar to FDR)

  genes <- rownames(rp_input$expr)

  results_df <- data.frame(
    gene_id = genes,
    # Up-regulation (class 1 > class 0)
    rp_up = rp_result$RPs[, 1],
    pfp_up = rp_result$pfp[, 1],
    pval_up = rp_result$pval[, 1],
    # Down-regulation (class 1 < class 0)
    rp_down = rp_result$RPs[, 2],
    pfp_down = rp_result$pfp[, 2],
    pval_down = rp_result$pval[, 2],
    stringsAsFactors = FALSE
  )

  # Combine into single logFC-like metric
  # Negative pfp_up = upregulated, positive pfp_down = downregulated
  results_df$direction <- ifelse(
    results_df$pfp_up < results_df$pfp_down, "up", "down"
  )
  results_df$pfp <- ifelse(
    results_df$direction == "up",
    results_df$pfp_up,
    results_df$pfp_down
  )
  results_df$pvalue <- ifelse(
    results_df$direction == "up",
    results_df$pval_up,
    results_df$pval_down
  )

  results_df <- results_df[order(results_df$pfp), ]

  cat("\nRankProd complete.\n")
  cat("Genes analyzed:", nrow(results_df), "\n")
  cat("Significant (PFP < 0.05):", sum(results_df$pfp < 0.05, na.rm = TRUE), "\n")
  cat("  Upregulated:", sum(results_df$pfp < 0.05 & results_df$direction == "up", na.rm = TRUE), "\n")
  cat("  Downregulated:", sum(results_df$pfp < 0.05 & results_df$direction == "down", na.rm = TRUE), "\n")

  return(results_df)
}


# ============================================================================
# Results Combination
# ============================================================================

#' Combine DExMA and RankProd results
#'
#' @param dexma_results Results from run_dexma_meta
#' @param rankprod_results Results from run_rankprod_meta
#' @param fdr_threshold FDR/PFP threshold
#' @return Combined results data frame
#' @export
combine_meta_results <- function(dexma_results,
                                  rankprod_results,
                                  fdr_threshold = 0.05) {

  # Merge results
  combined <- merge(
    dexma_results[, c("gene_id", "logFC", "fdr")],
    rankprod_results[, c("gene_id", "direction", "pfp")],
    by = "gene_id",
    all = TRUE,
    suffixes = c("_dexma", "_rankprod")
  )

  # Classify genes
  combined$sig_dexma <- !is.na(combined$fdr) & combined$fdr < fdr_threshold
  combined$sig_rankprod <- !is.na(combined$pfp) & combined$pfp < fdr_threshold

  combined$category <- "not_significant"
  combined$category[combined$sig_dexma & combined$sig_rankprod] <- "both"
  combined$category[combined$sig_dexma & !combined$sig_rankprod] <- "dexma_only"
  combined$category[!combined$sig_dexma & combined$sig_rankprod] <- "rankprod_only"

  # Direction agreement
  combined$direction_dexma <- ifelse(combined$logFC > 0, "up", "down")
  combined$direction_agree <- combined$direction_dexma == combined$direction

  cat("\n=== Combined Meta-Analysis Results ===\n\n")
  cat("Total genes:", nrow(combined), "\n")
  cat("Significant in both methods:", sum(combined$category == "both"), "\n")
  cat("DExMA only:", sum(combined$category == "dexma_only"), "\n")
  cat("RankProd only:", sum(combined$category == "rankprod_only"), "\n")
  cat("Direction agreement (in both):",
      sum(combined$direction_agree[combined$category == "both"], na.rm = TRUE),
      "/", sum(combined$category == "both"), "\n")

  return(combined)
}


#' Print summary of meta-analysis results
#'
#' @param results Meta-analysis results data frame
#' @param method Method name ("DExMA" or "RankProd")
#' @param fdr_col Column name for FDR/PFP values
#' @export
print_meta_summary <- function(results, method, fdr_col = "fdr") {
  cat("\n=======================================================\n")
  cat("           ", method, " RESULTS SUMMARY\n")
  cat("=======================================================\n\n")

  cat("Genes analyzed:", nrow(results), "\n\n")

  thresholds <- c(0.01, 0.05, 0.1)
  cat("Significant genes by threshold:\n")
  for (t in thresholds) {
    n_sig <- sum(results[[fdr_col]] < t, na.rm = TRUE)
    cat(sprintf("  FDR < %.2f: %d genes\n", t, n_sig))
  }
  cat("\n")
}


# ============================================================================
# Metafor Meta-Analysis (Alternative - no special system dependencies)
# ============================================================================

#' Run metafor-based effect-size meta-analysis
#'
#' Uses limma for per-study DE, then metafor for random-effects meta-analysis.
#' Works without DExMA/RankProd system dependencies.
#'
#' @param studies List of study objects from load_studies_for_meta
#' @param min_studies Minimum studies a gene must appear in
#' @return Data frame with meta-analysis results
#' @export
run_metafor_meta <- function(studies, min_studies = 4) {

  if (!requireNamespace("metafor", quietly = TRUE)) {
    stop("metafor package not installed. Run: install.packages('metafor')")
  }

  cat("\n=== Running Metafor Effect-Size Meta-Analysis ===\n\n")
  cat("Minimum studies per gene:", min_studies, "\n")
  cat("Studies:", length(studies), "\n\n")

  # Run per-study limma
  cat("Running per-study differential expression...\n")
  per_study_de <- list()

  for (study_id in names(studies)) {
    study <- studies[[study_id]]
    de <- run_limma_single_study(study$expr, study$groups)

    # Add sample sizes for effect size calculation
    de$n1 <- study$n_contrast
    de$n0 <- study$n_baseline
    de$study <- study_id

    per_study_de[[study_id]] <- de
    cat("  ", study_id, ":", nrow(de), "genes\n")
  }

  # Get genes present in enough studies
  all_genes <- unique(unlist(lapply(per_study_de, function(d) d$gene_id)))
  gene_counts <- table(unlist(lapply(per_study_de, function(d) d$gene_id)))
  genes_keep <- names(gene_counts[gene_counts >= min_studies])

  cat("\nGenes in >=", min_studies, "studies:", length(genes_keep), "\n")
  cat("Running meta-analysis per gene...\n\n")

  # Meta-analyze each gene
  results <- data.frame(
    gene_id = character(),
    logFC = numeric(),
    se = numeric(),
    pvalue = numeric(),
    n_studies = integer(),
    stringsAsFactors = FALSE
  )

  pb_interval <- max(1, length(genes_keep) %/% 20)

  for (i in seq_along(genes_keep)) {
    gene <- genes_keep[i]

    if (i %% pb_interval == 0) {
      cat(sprintf("  Progress: %d/%d (%.0f%%)\n",
                  i, length(genes_keep), 100 * i / length(genes_keep)))
    }

    # Collect effect sizes from each study
    yi <- c()  # Effect sizes (logFC)
    vi <- c()  # Variances

    for (study_id in names(per_study_de)) {
      de <- per_study_de[[study_id]]
      gene_row <- de[de$gene_id == gene, ]

      if (nrow(gene_row) == 1) {
        # Use logFC as effect size, SE^2 as variance
        yi <- c(yi, gene_row$logFC)
        # Variance = SE^2, and SE is available from limma via t-statistic
        # SE = logFC / t
        se <- abs(gene_row$logFC / gene_row$t)
        if (is.finite(se) && se > 0) {
          vi <- c(vi, se^2)
        } else {
          # Fallback: estimate variance from sample sizes
          n_total <- gene_row$n0 + gene_row$n1
          vi <- c(vi, 4 / n_total)  # Approximate variance
        }
      }
    }

    if (length(yi) >= min_studies) {
      # Random-effects meta-analysis
      tryCatch({
        meta_fit <- metafor::rma(yi = yi, vi = vi, method = "REML")

        results <- rbind(results, data.frame(
          gene_id = gene,
          logFC = meta_fit$b[1],
          se = meta_fit$se,
          pvalue = meta_fit$pval,
          n_studies = length(yi),
          stringsAsFactors = FALSE
        ))
      }, error = function(e) {
        # Skip genes where meta-analysis fails
      })
    }
  }

  # Adjust p-values
  results$fdr <- p.adjust(results$pvalue, method = "BH")
  results <- results[order(results$fdr), ]

  cat("\nMetafor meta-analysis complete.\n")
  cat("Genes analyzed:", nrow(results), "\n")
  cat("Significant (FDR < 0.05):", sum(results$fdr < 0.05, na.rm = TRUE), "\n")
  cat("  Upregulated:", sum(results$fdr < 0.05 & results$logFC > 0, na.rm = TRUE), "\n")
  cat("  Downregulated:", sum(results$fdr < 0.05 & results$logFC < 0, na.rm = TRUE), "\n")

  return(results)
}
