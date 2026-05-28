#!/usr/bin/env Rscript
# Bar plot of sample counts by gestational week for a config run,
# colored by exact vs estimated gestational age.
# Uses the phase2b phenodata TSV directly (not samples.csv).

source("scripts/visualization/ga_helpers.R")

plot_ga_exact_vs_estimated <- function(config_path, output_path = NULL,
                                       width = 1400, height = 800,
                                       res = 120) {
  config <- yaml::read_yaml(config_path)
  datasets <- config$files$datasets
  conditions <- tolower(config$sample_filter$condition)
  trimesters <- config$sample_filter$trimester

  pheno <- read.delim(config$paths$phenodata, stringsAsFactors = FALSE)

  col_map <- config$paths$column_map
  pheno$dataset <- pheno[[col_map$secondaryaccession]]
  pheno$sample_id <- pheno[[col_map$arraydatafile_exprscolumnnames]]
  pheno$trimester <- pheno[[col_map$Gestational.Age.Category]]
  pheno$condition <- pheno[[col_map$Diagnosis]]

  pheno <- pheno[pheno$dataset %in% unlist(datasets) &
                 tolower(pheno$condition) %in% conditions &
                 pheno$trimester %in% trimesters, ]

  exprs_file <- file.path(config$paths$output, "exprs_softimpute_combat_ref.tsv")
  exprs_cols <- colnames(read.delim(exprs_file, row.names = 1,
                                     check.names = FALSE, nrows = 1))
  pheno <- pheno[pheno$sample_id %in% exprs_cols, ]

  cat(sprintf("Samples after filtering: %d\n", nrow(pheno)))
  cat(sprintf("  Datasets: %s\n", paste(sort(unique(pheno$dataset)), collapse = ", ")))
  print(table(pheno$trimester))

  exact_week <- pheno$fetus_week
  exact_week[is.na(exact_week)] <- 0L

  range_weeks <- as.character(pheno$fetus_range_week)
  range_weeks[is.na(range_weeks) | range_weeks == "0"] <- ""

  df <- data.frame(
    sample_id = pheno$sample_id,
    dataset = pheno$dataset,
    trimester = pheno$trimester,
    exact_week = exact_week,
    range_weeks = range_weeks,
    rep_block = NA_character_,
    stringsAsFactors = FALSE
  )

  df <- fill_range_from_category(df)
  df <- estimate_weeks(df)

  exact <- df[df$exact_week > 0, ]
  est <- df[df$exact_week == 0, ]

  week_min <- min(df$estimated_week, na.rm = TRUE)
  week_max <- max(df$estimated_week, na.rm = TRUE)
  all_weeks <- seq(week_min, week_max)
  wk_labels <- as.character(all_weeks)

  exact_counts <- table(factor(exact$estimated_week, levels = all_weeks))
  est_counts <- table(factor(est$estimated_week, levels = all_weeks))

  count_mat <- rbind(Exact = as.integer(exact_counts),
                     Estimated = as.integer(est_counts))
  colnames(count_mat) <- wk_labels

  keep <- colSums(count_mat) > 0
  count_mat <- count_mat[, keep, drop = FALSE]

  col_exact <- "#56B4E9"
  col_est <- "#E69F00"

  n_exact <- nrow(exact)
  n_est <- nrow(est)
  subtitle <- sprintf("%d exact + %d estimated weeks (%d samples total)",
                      n_exact, n_est, n_exact + n_est)

  do_plot <- function() {
    par(mar = c(5, 5, 4, 2))
    ymax <- max(colSums(count_mat)) * 1.15

    bp <- barplot(count_mat,
                  col = c(col_exact, col_est),
                  border = NA,
                  xlab = "Gestational week",
                  ylab = "Number of samples",
                  main = "Sample distribution by gestational week",
                  cex.main = 1.3, cex.lab = 1.1, cex.axis = 1.0,
                  las = 1, ylim = c(0, ymax))
    mtext(subtitle, side = 3, line = 0.3, cex = 0.8, col = "grey40")

    t2_col <- which(colnames(count_mat) == "13")
    if (length(t2_col) == 1 && length(bp) > 1) {
      x_border <- bp[t2_col] - diff(bp[1:2]) / 2
      abline(v = x_border, lty = 2, col = "grey50", lwd = 1.5)
      text(x_border, ymax * 0.95, "T1 | T2",
           cex = 0.8, col = "grey50")
    }

    legend("topright",
           legend = c(sprintf("Exact week (n=%d)", n_exact),
                      sprintf("Estimated (n=%d)", n_est)),
           fill = c(col_exact, col_est),
           border = NA, cex = 0.9, bty = "n")
  }

  if (!is.null(output_path)) {
    dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
    png(output_path, width = width, height = height, res = res)
    do_plot()
    dev.off()
    cat("Saved:", output_path, "\n")
  } else {
    do_plot()
  }

  invisible(list(exact = n_exact, estimated = n_est, count_mat = count_mat))
}

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  config_path <- if (length(args) >= 1) args[1] else
    paste0("scripts/integrative_analysis/phase2b_direct_merge/",
           "config_yehor_sashko/",
           "config_phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko.yaml")
  output_path <- if (length(args) >= 2) args[2] else
    "articles/yehor_conference_2026/data/fig_ga_exact_vs_estimated_6ds.png"
  plot_ga_exact_vs_estimated(config_path, output_path)
}
