#!/usr/bin/env Rscript
#' Plot NA staircase diagrams for all non-archived phase2b runs
#' on a single figure with shared axis scales — one PNG per FDR method.

suppressPackageStartupMessages(library(yaml))

script_dir <- if (exists("script_dir")) {
  script_dir
} else {
  "scripts/integrative_analysis/phase2b_direct_merge"
}
source(file.path(script_dir, "plot_na_staircase.R"))

setwd("/home/shivers/a/r/igea-r/scripts/integrative-gene-expression-analysis")

config_dir <- "scripts/integrative_analysis/phase2b_direct_merge"
yamls <- sort(list.files(config_dir, pattern = "^config_phase2b_.*\\.yaml$",
                         full.names = TRUE, recursive = TRUE))

runs <- list()
for (yf in yamls) {
  cfg <- yaml::read_yaml(yf)
  out_dir <- cfg$paths$output
  if (!dir.exists(out_dir)) next
  label <- sub("^config_phase2b_", "", sub("\\.yaml$", "", basename(yf)))
  merged <- build_merged_from_config(yf)
  runs[[label]] <- list(merged = merged, label = label)
}

if (length(runs) == 0) stop("No non-archived runs found")
cat(sprintf("Found %d runs with existing output\n", length(runs)))

all_methods <- unique(unlist(lapply(runs, function(r) names(r$merged$gene_fdr_list))))
if (length(all_methods) == 0) all_methods <- "none"

n <- length(runs)
ncols <- min(n, 3)
nrows <- ceiling(n / ncols)

for (m in all_methods) {
  staircases <- lapply(runs, function(r) {
    gf <- r$merged$gene_fdr_list[[m]]
    prepare_staircase(r$merged$matrix, gene_fdr = gf)
  })

  max_genes   <- max(sapply(staircases, function(s) s$nr_orig))
  max_samples <- max(sapply(staircases, function(s) s$nc))

  for (i in seq_along(staircases)) {
    s <- staircases[[i]]
    deficit <- max_genes - s$nr_orig
    if (deficit > 0) {
      pad_na <- matrix(TRUE, deficit, s$nc)
      s$sorted_na  <- rbind(pad_na, s$sorted_na)
      s$excluded   <- c(rep(FALSE, deficit), s$excluded)
      s$gene_names <- c(rep(NA_character_, deficit), s$gene_names)
      s$actual_genes <- s$nr_orig
      s$nr_orig      <- max_genes
    } else {
      s$actual_genes <- s$nr_orig
    }
    staircases[[i]] <- s
  }

  out_png <- sprintf("output/na_staircase_all_runs_%s.png", m)
  png(out_png, width = 900 * ncols, height = 700 * nrows, res = 150)
  par(mfrow = c(nrows, ncols), mar = c(4, 5, 3, 1), oma = c(0, 0, 2, 0))

  for (i in seq_along(runs)) {
    render_staircase(
      staircases[[i]],
      title      = runs[[i]]$label,
      sample_ds  = runs[[i]]$merged$sample_ds,
      xlim       = c(0, max_samples),
      cex_main   = 0.9,
      cex_legend = 0.45
    )
  }

  remaining <- ncols * nrows - n
  if (remaining > 0) {
    for (i in seq_len(remaining)) plot.new()
  }

  mtext(sprintf("NA Missingness \u2014 All Phase 2B Runs [FDR: %s]", m),
        outer = TRUE, cex = 1.0, line = 0.5)
  dev.off()
  cat("Wrote:", out_png, "\n")
}
