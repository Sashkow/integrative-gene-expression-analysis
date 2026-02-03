#!/usr/bin/env Rscript

#' Run Phase 2: Meta-Analysis
#'
#' Usage:
#'   Rscript run_phase2.R [options]
#'
#' Options:
#'   --config=PATH         Path to config YAML (default: config_phase2.yaml)
#'   --method=NAME         Run only specific method: dexma, rankprod, or both
#'   --output_dir=PATH     Override output directory
#'   --no_archive          Don't archive previous results
#'   --help                Show this help message

library(yaml)

# Null coalescing operator
`%||%` <- function(a, b) if (is.null(a)) b else a

# Source shared utilities
utils_path <- "scripts/utils/logging_utils.R"
if (!file.exists(utils_path)) {
  utils_path <- file.path("..", "..", "utils", "logging_utils.R")
}
source(utils_path)

# Get script directory
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("--file=", "", file_arg))))
  }
  return("scripts/integrative_analysis/phase2_meta_analysis")
}

script_dir <- get_script_dir()

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Defaults
config_path <- file.path(script_dir, "config_phase2.yaml")
override_method <- NULL
override_output <- NULL
no_archive <- FALSE

# Parse arguments
for (arg in args) {
  if (arg == "--help" || arg == "-h") {
    cat("
Phase 2: Meta-Analysis

Usage:
  Rscript run_phase2.R [options]

Options:
  --config=PATH         Path to config YAML (default: config_phase2.yaml)
  --method=NAME         Run only specific method: dexma, rankprod, or both
  --output_dir=PATH     Override output directory
  --no_archive          Don't archive previous results
  --help                Show this help message

Examples:
  Rscript run_phase2.R
  Rscript run_phase2.R --method=dexma
  Rscript run_phase2.R --method=rankprod
")
    quit(status = 0)
  } else if (grepl("^--config=", arg)) {
    config_path <- sub("^--config=", "", arg)
  } else if (grepl("^--method=", arg)) {
    override_method <- sub("^--method=", "", arg)
  } else if (grepl("^--output_dir=", arg)) {
    override_output <- sub("^--output_dir=", "", arg)
  } else if (arg == "--no_archive") {
    no_archive <- TRUE
  }
}

# Load config
if (!file.exists(config_path)) {
  stop("Config file not found: ", config_path)
}
config <- yaml::read_yaml(config_path)

# Apply overrides
if (!is.null(override_output)) config$paths$output <- override_output
if (!is.null(override_method)) {
  if (override_method == "dexma") {
    config$methods$dexma$enabled <- TRUE
    config$methods$rankprod$enabled <- FALSE
  } else if (override_method == "rankprod") {
    config$methods$dexma$enabled <- FALSE
    config$methods$rankprod$enabled <- TRUE
  }
}

# Extract config values
mapped_path <- config$paths$mapped_data
phenodata_path <- config$paths$phenodata
output_dir <- config$paths$output
datasets <- config$files$datasets

# Phenotype settings
group_column <- config$phenotype$group_column
baseline <- config$phenotype$baseline
contrast <- config$phenotype$contrast

# Method settings
run_metafor <- config$methods$metafor$enabled %||% TRUE
run_dexma <- config$methods$dexma$enabled %||% FALSE
run_rankprod <- config$methods$rankprod$enabled %||% FALSE

# Thresholds
fdr_threshold <- config$thresholds$fdr %||% 0.05
logfc_threshold <- config$thresholds$logfc %||% 1.0

cat("\n")
cat("============================================================\n")
cat("  Phase 2: Meta-Analysis\n")
cat("============================================================\n\n")

cat("Configuration:\n")
cat("  Config file:  ", config_path, "\n")
cat("  Mapped path:  ", mapped_path, "\n")
cat("  Phenodata:    ", phenodata_path, "\n")
cat("  Output dir:   ", output_dir, "\n")
cat("  Comparison:   ", contrast, "vs", baseline, "\n")
if (!is.null(datasets)) {
  cat("  Datasets:     ", paste(datasets, collapse = ", "), "\n")
}
cat("\n")
cat("Methods to run:\n")
cat("  Metafor:      ", run_metafor, "\n")
cat("  DExMA:        ", run_dexma, "\n")
cat("  RankProd:     ", run_rankprod, "\n")
cat("\n")

# Source the module
source(file.path(script_dir, "meta_analysis.R"))

# Archive previous results
if (!no_archive) {
  archive_previous_results(output_dir)
} else {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
}

# Setup logging
log_file <- NULL
if (config$logging$save_log %||% TRUE) {
  log_file <- setup_logging(output_dir, prefix = "phase2_log")
}

# Load studies
cat("Loading studies...\n\n")
studies <- load_studies_for_meta(
  mapped_path = mapped_path,
  phenodata_path = phenodata_path,
  datasets = datasets,
  group_column = group_column,
  baseline = baseline,
  contrast = contrast
)

cat("\nLoaded", length(studies), "studies\n")

# Filter to protein-coding if configured
if (config$gene_filter$protein_coding_only %||% FALSE) {
  cat("\nFiltering to protein-coding genes...\n")
  studies <- filter_studies_protein_coding(studies)
}

# Run per-study DE if configured
per_study_de <- NULL
if (config$output$save_per_study_de %||% FALSE) {
  cat("\n=== Running Per-Study Differential Expression ===\n\n")
  per_study_de <- run_per_study_de(studies)

  # Save per-study results
  for (study_id in names(per_study_de)) {
    de_file <- file.path(output_dir, paste0(study_id, "_de.tsv"))
    write.table(per_study_de[[study_id]], de_file, sep = "\t",
                row.names = FALSE, quote = FALSE)
  }
  cat("\nSaved per-study DE results\n")
}

# Run Metafor (default method)
metafor_results <- NULL
if (run_metafor) {
  metafor_results <- tryCatch({
    run_metafor_meta(
      studies,
      min_studies = config$gene_filter$min_studies %||% 4
    )
  }, error = function(e) {
    cat("ERROR running metafor:", conditionMessage(e), "\n")
    NULL
  })

  if (!is.null(metafor_results)) {
    metafor_file <- file.path(output_dir, "metafor_results.tsv")
    write.table(metafor_results, metafor_file, sep = "\t",
                row.names = FALSE, quote = FALSE)
    cat("Saved metafor results to:", metafor_file, "\n")

    # Save significant genes
    metafor_sig <- metafor_results[metafor_results$fdr < fdr_threshold, ]
    if (nrow(metafor_sig) > 0) {
      metafor_sig_file <- file.path(output_dir, "metafor_significant.tsv")
      write.table(metafor_sig, metafor_sig_file, sep = "\t",
                  row.names = FALSE, quote = FALSE)
      cat("Saved", nrow(metafor_sig), "significant genes to:",
          metafor_sig_file, "\n")
    }

    print_meta_summary(metafor_results, "Metafor", "fdr")
  }
}

# Run DExMA (requires balanced studies with samples in both groups)
dexma_results <- NULL
if (run_dexma) {
  balanced_studies <- filter_balanced_studies(studies, min_per_group = 2)

  if (length(balanced_studies) < 2) {
    cat("WARNING: DExMA requires at least 2 balanced studies. Found:",
        length(balanced_studies), "\n")
    cat("Skipping DExMA analysis.\n")
  } else {
    dexma_results <- tryCatch({
      run_dexma_meta(
        balanced_studies,
      effect_size = config$methods$dexma$effect_size %||% "SMD",
      missAllow = config$methods$dexma$missAllow %||% 0.3,
      impute = config$methods$dexma$imputation$enabled %||% TRUE
    )
  }, error = function(e) {
    cat("ERROR running DExMA:", conditionMessage(e), "\n")
    NULL
  })

  if (!is.null(dexma_results) && (config$output$save_dexma_results %||% TRUE)) {
    dexma_file <- file.path(output_dir, "dexma_results.tsv")
    write.table(dexma_results, dexma_file, sep = "\t",
                row.names = FALSE, quote = FALSE)
    cat("Saved DExMA results to:", dexma_file, "\n")

    # Save significant genes
    dexma_sig <- dexma_results[dexma_results$fdr < fdr_threshold, ]
    if (nrow(dexma_sig) > 0) {
      dexma_sig_file <- file.path(output_dir, "dexma_significant.tsv")
      write.table(dexma_sig, dexma_sig_file, sep = "\t",
                  row.names = FALSE, quote = FALSE)
      cat("Saved", nrow(dexma_sig), "significant genes to:", dexma_sig_file, "\n")
    }

    print_meta_summary(dexma_results, "DExMA", "fdr")
    }
  }
}

# Run RankProd (requires balanced studies with samples in both groups)
rankprod_results <- NULL
if (run_rankprod) {
  balanced_studies_rp <- filter_balanced_studies(studies, min_per_group = 2)

  if (length(balanced_studies_rp) < 2) {
    cat("WARNING: RankProd requires at least 2 balanced studies. Found:",
        length(balanced_studies_rp), "\n")
    cat("Skipping RankProd analysis.\n")
  } else {
    rankprod_results <- tryCatch({
      run_rankprod_meta(
        balanced_studies_rp,
        min_studies = min(config$gene_filter$min_studies %||% 4,
                         length(balanced_studies_rp)),
        num_perm = config$methods$rankprod$num_perm %||% 1000,
        logged = config$methods$rankprod$logged %||% TRUE
      )
    }, error = function(e) {
      cat("ERROR running RankProd:", conditionMessage(e), "\n")
      NULL
    })

    if (!is.null(rankprod_results) &&
        (config$output$save_rankprod_results %||% TRUE)) {
      rankprod_file <- file.path(output_dir, "rankprod_results.tsv")
      write.table(rankprod_results, rankprod_file, sep = "\t",
                  row.names = FALSE, quote = FALSE)
      cat("Saved RankProd results to:", rankprod_file, "\n")

      # Save significant genes
      rankprod_sig <- rankprod_results[rankprod_results$pfp < fdr_threshold, ]
      if (nrow(rankprod_sig) > 0) {
        rankprod_sig_file <- file.path(output_dir, "rankprod_significant.tsv")
        write.table(rankprod_sig, rankprod_sig_file, sep = "\t",
                    row.names = FALSE, quote = FALSE)
        cat("Saved", nrow(rankprod_sig), "significant genes to:",
            rankprod_sig_file, "\n")
      }

      print_meta_summary(rankprod_results, "RankProd", "pfp")
    }
  }
}

# Combine results
if (!is.null(dexma_results) && !is.null(rankprod_results) &&
    (config$output$save_combined_results %||% TRUE)) {

  combined <- combine_meta_results(dexma_results, rankprod_results, fdr_threshold)

  combined_file <- file.path(output_dir, "combined_results.tsv")
  write.table(combined, combined_file, sep = "\t",
              row.names = FALSE, quote = FALSE)
  cat("\nSaved combined results to:", combined_file, "\n")

  # Save high-confidence genes (significant in both)
  both_sig <- combined[combined$category == "both", ]
  if (nrow(both_sig) > 0) {
    both_file <- file.path(output_dir, "high_confidence_degs.tsv")
    write.table(both_sig, both_file, sep = "\t",
                row.names = FALSE, quote = FALSE)
    cat("Saved", nrow(both_sig), "high-confidence DEGs to:", both_file, "\n")
  }
}

# Save summary
summary_file <- file.path(output_dir, "summary.txt")
summary_conn <- file(summary_file, "w")
writeLines(c(
  "Phase 2: Meta-Analysis Summary",
  "==============================",
  "",
  paste("Date:", Sys.time()),
  paste("Config:", config_path),
  "",
  paste("Comparison:", contrast, "vs", baseline),
  paste("Studies analyzed:", length(studies)),
  "",
  "Results:",
  if (!is.null(metafor_results)) {
    c(
      paste("  Metafor genes:", nrow(metafor_results)),
      paste("  Metafor significant (FDR <", fdr_threshold, "):",
            sum(metafor_results$fdr < fdr_threshold, na.rm = TRUE))
    )
  } else {
    "  Metafor: not run"
  },
  if (!is.null(dexma_results)) {
    c(
      paste("  DExMA genes:", nrow(dexma_results)),
      paste("  DExMA significant (FDR <", fdr_threshold, "):",
            sum(dexma_results$fdr < fdr_threshold, na.rm = TRUE))
    )
  } else {
    "  DExMA: not run"
  },
  if (!is.null(rankprod_results)) {
    c(
      paste("  RankProd genes:", nrow(rankprod_results)),
      paste("  RankProd significant (PFP <", fdr_threshold, "):",
            sum(rankprod_results$pfp < fdr_threshold, na.rm = TRUE))
    )
  } else {
    "  RankProd: not run"
  }
), summary_conn)
close(summary_conn)
cat("\nSaved summary to:", summary_file, "\n")

# Close logging
if (!is.null(log_file)) {
  close_logging(log_file)
}

cat("\n============================================================\n")
cat("  Phase 2 Complete!\n")
cat("============================================================\n\n")

cat("Output files in", output_dir, ":\n")
for (f in list.files(output_dir, pattern = "\\.(tsv|txt)$")) {
  cat("  ", f, "\n")
}
cat("\n")
