#!/usr/bin/env Rscript
# Stacked barplot of sample counts by gestational week, colored by dataset.
# Samples with exact weeks are solid; range-only samples are striped.

source("scripts/visualization/ga_helpers.R")

plot_gestational_age_by_dataset <- function(config_path, output_path = NULL,
                                            width = 1400, height = 800,
                                            res = 120) {
  config <- yaml::read_yaml(config_path)
  datasets <- config$files$datasets
  conditions <- tolower(config$sample_filter$condition)
  trimesters <- config$sample_filter$trimester

  samples_csv <- read.csv("data/phenodata/samples.csv",
                          stringsAsFactors = FALSE)

  used_samples <- character(0)
  for (ds in datasets) {
    fname <- config$files$file_map[[ds]]
    if (is.null(fname))
      fname <- paste0(ds, config$files$suffix, ".tsv")
    fpath <- file.path(config$paths$mapped_data, fname)
    mat <- read.delim(fpath, row.names = 1,
                      check.names = FALSE, nrows = 1)
    ds_samples <- colnames(mat)
    idx <- which(
      samples_csv$secondaryaccession == ds &
      samples_csv$arraydatafile_exprscolumnnames %in% ds_samples &
      tolower(samples_csv$Diagnosis) %in% conditions &
      tolower(samples_csv$Gestational.Age.Category) %in%
        tolower(trimesters)
    )
    used_samples <- c(used_samples,
                      samples_csv$arraydatafile_exprscolumnnames[idx])
  }

  pheno <- samples_csv[
    samples_csv$arraydatafile_exprscolumnnames %in% used_samples, ]

  df <- parse_ga_rows(pheno)
  df <- collapse_rep_blocks(df)
  df <- fill_range_from_category(df)
  df <- estimate_weeks(df)

  # --- Build count matrices ---
  all_ds <- sort(unique(df$dataset))
  exact <- df[df$exact_week > 0, ]
  est <- df[df$exact_week == 0, ]

  week_min <- min(df$estimated_week, na.rm = TRUE)
  week_max <- max(df$estimated_week, na.rm = TRUE)
  all_weeks <- seq(week_min, week_max)

  exact_mat <- make_count_matrix(exact$dataset, exact$exact_week,
                                 all_ds, all_weeks)
  est_mat <- make_count_matrix(est$dataset, est$estimated_week,
                               all_ds, all_weeks)

  keep_weeks <- colSums(exact_mat) + colSums(est_mat) > 0
  exact_mat <- exact_mat[, keep_weeks, drop = FALSE]
  est_mat <- est_mat[, keep_weeks, drop = FALSE]

  oi_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                  "#0072B2", "#D55E00", "#CC79A7", "#999999")
  ds_colors <- setNames(oi_palette[seq_along(all_ds)], all_ds)

  n_exact <- nrow(exact)
  n_est <- nrow(est)
  subtitle <- sprintf(
    "%d exact + %d estimated (evenly distributed across known range)",
    n_exact, n_est)

  do_plot <- function() {
    par(mar = c(5, 5, 4, 2))
    combined_totals <- colSums(exact_mat) + colSums(est_mat)
    ymax <- max(combined_totals) * 1.15

    bp <- barplot(exact_mat, col = ds_colors[all_ds],
                  border = NA,
                  xlab = "Gestational week",
                  ylab = "Number of samples",
                  main = "Sample distribution by gestational week",
                  cex.main = 1.3, cex.lab = 1.1, cex.axis = 1.0,
                  las = 1, ylim = c(0, ymax))
    mtext(subtitle, side = 3, line = 0.3, cex = 0.8, col = "grey40")

    bar_width <- if (length(bp) > 1) diff(bp[1:2]) * 0.8 else bp[1] * 0.8

    for (wi in seq_along(bp)) {
      base <- colSums(exact_mat)[wi]
      for (di in seq_along(all_ds)) {
        n_here <- est_mat[di, wi]
        if (n_here == 0) next
        x0 <- bp[wi] - bar_width / 2
        x1 <- bp[wi] + bar_width / 2
        y0 <- base; y1 <- base + n_here
        rect(x0, y0, x1, y1, col = ds_colors[all_ds[di]],
             border = NA)
        stripe_step <- 0.3
        clip(x0, x1, y0, y1)
        span <- (x1 - x0) + (y1 - y0)
        for (s in seq(-span, span, by = stripe_step))
          lines(c(x0, x0 + span), c(y0 + s, y0 + s + span),
                col = "white", lwd = 1.2)
        do.call("clip", as.list(par("usr")))
        base <- y1
      }
    }

    t2_col <- which(colnames(exact_mat) == "13")
    if (length(t2_col) == 1) {
      x_border <- bp[t2_col] - diff(bp[1:2]) / 2
      abline(v = x_border, lty = 2, col = "grey50", lwd = 1.5)
      text(x_border, ymax * 0.95, "T1 | T2",
           cex = 0.8, col = "grey50")
    }

    exact_n <- rowSums(exact_mat)
    est_n <- rowSums(est_mat)
    has_samples <- (exact_n + est_n) > 0
    show_ds <- all_ds[has_samples]

    legend("topright",
           legend = c(sprintf("%s (n=%d)", show_ds,
                              exact_n[show_ds] + est_n[show_ds]),
                      "", "solid = exact week",
                      "striped = estimated"),
           fill = c(ds_colors[show_ds], NA, "grey60", "grey60"),
           density = c(rep(NA, length(show_ds)), NA, NA, 20),
           angle = c(rep(0, length(show_ds)), 0, 0, 45),
           border = c(rep("white", length(show_ds)),
                      NA, "grey40", "grey40"),
           cex = 0.65, bty = "n", title = "Dataset")
  }

  if (!is.null(output_path)) {
    png(output_path, width = width, height = height, res = res)
    do_plot()
    dev.off()
    cat("Saved:", output_path, "\n")
  } else {
    do_plot()
  }

  invisible(list(exact = exact_mat, estimated = est_mat))
}

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  config_path <- if (length(args) >= 1) args[1] else
    paste0("scripts/integrative_analysis/phase2b_direct_merge/",
           "config_yehor_sashko/",
           "config_phase2b_1_2_yehor_7ds_no_37653_enriched_sashko.yaml")
  output_path <- if (length(args) >= 2) args[2] else
    "articles/imputation_article/figures/fig_gestational_age_by_dataset.png"
  plot_gestational_age_by_dataset(config_path, output_path)
}
