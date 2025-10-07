#' Differential Expression Analysis Module
#'
#' Functions for differential expression analysis using limma
#'
#' @author Expression Integration Pipeline
#' @date 2025

#' Perform differential expression analysis
#'
#' @param exprs Expression matrix
#' @param pdata Phenodata
#' @param design_formula Formula for design matrix
#' @param contrasts Character vector of contrasts
#' @param method Method for decideTests (default: "global")
#' @return List with results and fitted model
#' @export
differential_expression <- function(exprs, pdata,
                                    design_formula,
                                    contrasts,
                                    method = "global") {

  library(limma)

  # Create design matrix
  design <- model.matrix(design_formula, data = pdata)
  original_colnames <- colnames(design)
  message("Design matrix: ", ncol(design), " coefficients")
  message("  Original columns: ", paste(original_colnames, collapse = ", "))

  # Simplify column names for makeContrasts (it requires syntactically valid names)
  # Create mapping from simplified names back to originals
  simple_names <- paste0("C", seq_len(ncol(design)))
  colnames(design) <- simple_names
  name_mapping <- setNames(simple_names, original_colnames)

  message("  Simplified to: ", paste(simple_names, collapse = ", "))

  # Fit linear model
  fit <- lmFit(exprs, design)

  # Convert contrasts to use simplified names
  contrasts_fixed <- contrasts
  for (orig_name in names(name_mapping)) {
    contrasts_fixed <- gsub(orig_name, name_mapping[orig_name], contrasts_fixed, fixed = TRUE)
  }

  message("\nContrast translation:")
  for (i in seq_along(contrasts)) {
    message("  Original: ", contrasts[i])
    message("  Simplified: ", contrasts_fixed[i])
  }

  # Create contrast matrix
  contrast_matrix <- makeContrasts(contrasts = contrasts_fixed, levels = design)
  colnames(contrast_matrix) <- contrasts  # Use original names for output
  message("\nContrast matrix created successfully")

  # Fit contrasts
  fit2 <- contrasts.fit(fit, contrast_matrix)
  efit <- eBayes(fit2)

  # Decide tests
  results <- decideTests(efit, method = method)
  message("\nDifferentially expressed genes (method: ", method, "):")
  print(summary(results))

  return(list(
    efit = efit,
    results = results,
    contrast_matrix = contrast_matrix
  ))
}


#' Extract and annotate differential expression results
#'
#' @param efit Fitted eBayes object
#' @param results decideTests results
#' @param exprs Expression matrix
#' @param pdata Phenodata
#' @param species Species for annotation (default: "human")
#' @return Annotated differential expression table
#' @export
extract_difexp_results <- function(efit, results, exprs, pdata,
                                    species = "human") {

  library(limma)
  library(org.Hs.eg.db)
  library(AnnotationDbi)

  # Get top table
  tt <- topTable(efit, number = nrow(exprs), sort.by = "none")

  # Ensure order matches
  if (!all(rownames(results) == rownames(tt))) {
    stop("Row order mismatch between results and topTable")
  }

  # Combine results
  tt_results <- cbind(tt, as.data.frame(results))

  # Annotate with gene symbols and names
  if (species == "human") {
    annotation <- AnnotationDbi::select(
      org.Hs.eg.db,
      keys = rownames(tt_results),
      columns = c("SYMBOL", "GENENAME"),
      keytype = "ENTREZID"
    )

    # Remove duplicates and NAs
    annotation <- annotation[!duplicated(annotation$ENTREZID), ]
    annotation <- annotation[!is.na(annotation$SYMBOL), ]

    # Merge with results
    merged <- merge(annotation, tt_results,
                    by.x = "ENTREZID", by.y = "row.names",
                    all.y = FALSE)

  } else {
    merged <- tt_results
    merged$ENTREZID <- rownames(tt_results)
  }

  # Add average expression by group
  merged <- add_group_averages(merged, exprs, pdata)

  message("Differential expression results: ", nrow(merged), " genes annotated")

  return(merged)
}


#' Add group average expression values
#'
#' @param difexp Differential expression results
#' @param exprs Expression matrix
#' @param pdata Phenodata
#' @param group_col Column defining groups (default: "trim_term")
#' @return difexp with added average columns
#' @export
add_group_averages <- function(difexp, exprs, pdata,
                                group_col = "trim_term") {

  if (!group_col %in% colnames(pdata)) {
    warning("Group column '", group_col, "' not found")
    return(difexp)
  }

  # Get ENTREZID for indexing (expression matrix is indexed by ENTREZID)
  if (!"ENTREZID" %in% colnames(difexp)) {
    warning("No ENTREZID column found in differential expression results")
    return(difexp)
  }

  # Calculate averages for each group
  groups <- unique(as.character(pdata[[group_col]]))
  groups <- groups[!is.na(groups) & groups != ""]

  for (group in groups) {
    group_samples <- pdata[pdata[[group_col]] == group, ]$arraydatafile_exprscolumnnames
    group_samples <- make.names(group_samples)
    group_samples <- intersect(group_samples, colnames(exprs))

    if (length(group_samples) > 0) {
      col_name <- paste(group, "Average", sep = "_")

      # Use ENTREZID to index expression matrix
      # Only calculate for genes that exist in expression matrix
      gene_ids_in_exprs <- intersect(difexp$ENTREZID, rownames(exprs))

      # Initialize column with NA
      difexp[[col_name]] <- NA

      # Calculate averages for genes that exist in expression matrix
      idx <- match(gene_ids_in_exprs, difexp$ENTREZID)
      if (length(group_samples) == 1) {
        difexp[[col_name]][idx] <- exprs[gene_ids_in_exprs, group_samples]
      } else {
        difexp[[col_name]][idx] <- rowMeans(exprs[gene_ids_in_exprs, group_samples, drop = FALSE])
      }
      message("  Added column: ", col_name)
    }
  }

  return(difexp)
}


#' Filter differential expression results
#'
#' @param difexp Differential expression results
#' @param logfc_threshold Log fold-change threshold (default: 1)
#' @param pvalue_threshold Adjusted p-value threshold (default: 0.05)
#' @return Filtered differential expression results
#' @export
filter_difexp <- function(difexp, logfc_threshold = 1, pvalue_threshold = 0.05) {

  if ("logFC" %in% colnames(difexp) && "adj.P.Val" %in% colnames(difexp)) {
    filtered <- difexp[abs(difexp$logFC) >= logfc_threshold &
                       difexp$adj.P.Val < pvalue_threshold, ]
  } else if ("logFC" %in% colnames(difexp)) {
    filtered <- difexp[abs(difexp$logFC) >= logfc_threshold, ]
  } else {
    warning("No logFC column found")
    return(difexp)
  }

  message("Filtered: ", nrow(filtered), " genes (|logFC| >= ", logfc_threshold,
          ", adj.P.Val < ", pvalue_threshold, ")")

  return(filtered)
}
