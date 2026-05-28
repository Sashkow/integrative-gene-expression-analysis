#!/usr/bin/env Rscript
# Stacked barplot of sample counts by gestational week for a given tissue type.
# Config-driven via YAML: tissue types, condition, output location.
# Samples with exact weeks are solid; range-only samples are striped.
#
# Usage: Rscript plot_gestational_age_by_tissue.R <config.yaml>

library(openxlsx)
source("scripts/visualization/ga_helpers.R")

plot_gestational_age_by_tissue <- function(config_path,
                                           width = 2000, height = 1000,
                                           res = 120) {

  config <- yaml::read_yaml(config_path)
  tissues <- config$tissues
  condition <- tolower(config$condition)
  tissue_label <- config$tissue_label
  output_dir <- config$output_dir
  file_prefix <- config$file_prefix
  require_exact_for_term_only <- isTRUE(config$require_exact_for_term_only)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  samples_csv <- read.csv("data/phenodata/samples.csv", stringsAsFactors = FALSE)

  pheno <- samples_csv[tolower(samples_csv$Diagnosis) == condition &
                        samples_csv$Biological.Specimen %in% tissues &
                        !(samples_csv$X.excluded. %in% TRUE) &
                        !(samples_csv$secondaryaccession %in% c("", "_", NA)), ]

  if (nrow(pheno) == 0) {
    cat(sprintf("No %s %s samples found.\n", condition, tissue_label))
    return(invisible(NULL))
  }

  cat(sprintf("Total %s %s samples: %d\n", condition, tissue_label, nrow(pheno)))
  cat("Datasets:", paste(sort(unique(pheno$secondaryaccession)), collapse = ", "), "\n\n")

  df <- parse_ga_rows(pheno)
  df$condition <- condition
  df <- collapse_rep_blocks(df)

  if (require_exact_for_term_only) {
    term_categories <- c("Term", "Third trimester")
    drop_ds <- character(0)
    for (ds in unique(df$dataset)) {
      ds_rows <- df[df$dataset == ds, ]
      all_term <- all(ds_rows$trimester %in% term_categories)
      if (all_term) {
        has_exact <- any(ds_rows$exact_week > 0)
        has_range <- any(!is.na(ds_rows$range_weeks) & ds_rows$range_weeks != "")
        if (!has_exact && !has_range) {
          drop_ds <- c(drop_ds, ds)
        }
      }
    }
    if (length(drop_ds) > 0) {
      cat(sprintf("Excluded %d term-only datasets without exact/range GA: %s\n",
                  length(drop_ds), paste(drop_ds, collapse = ", ")))
      df <- df[!df$dataset %in% drop_ds, ]
    }
  }

  df <- fill_range_from_category(df)
  df <- estimate_weeks(df)

  if (nrow(df) == 0) {
    cat("No samples with estimable gestational age.\n")
    return(invisible(NULL))
  }

  # --- Count matrices ---
  all_ds <- sort(unique(df$dataset))
  exact <- df[df$exact_week > 0, ]
  est <- df[df$exact_week == 0, ]

  week_min <- min(df$estimated_week, na.rm = TRUE)
  week_max <- max(df$estimated_week, na.rm = TRUE)
  all_weeks <- seq(week_min, week_max)

  exact_mat <- make_count_matrix(exact$dataset, exact$exact_week, all_ds, all_weeks)
  est_mat <- make_count_matrix(est$dataset, est$estimated_week, all_ds, all_weeks)

  ds_colors <- setNames(GA_BASE_PALETTE[seq_along(all_ds)], all_ds)

  n_exact <- nrow(exact)
  n_est <- nrow(est)
  subtitle <- sprintf("%d exact + %d estimated from %d datasets (%s, %s)",
                      n_exact, n_est, length(all_ds), condition, tissue_label)

  # --- Plot ---
  plot_path <- file.path(output_dir, paste0(file_prefix, ".png"))
  png(plot_path, width = width, height = height, res = res)

  par(mar = c(5, 5, 4, 2))
  combined_totals <- colSums(exact_mat) + colSums(est_mat)
  ymax <- max(combined_totals) * 1.15

  bp <- barplot(exact_mat, col = ds_colors[all_ds],
                border = NA,
                xlab = "Gestational week", ylab = "Number of samples",
                main = sprintf("Sample distribution by gestational week (%s)", tissue_label),
                cex.main = 1.1, cex.lab = 1.0, cex.axis = 0.8,
                las = 2, ylim = c(0, ymax))
  mtext(subtitle, side = 3, line = 0.3, cex = 0.7, col = "grey40")

  bar_width <- if (length(bp) > 1) diff(bp[1:2]) * 0.8 else bp[1] * 0.8
  for (wi in seq_along(bp)) {
    base <- colSums(exact_mat)[wi]
    for (di in seq_along(all_ds)) {
      n_here <- est_mat[di, wi]
      if (n_here == 0) next
      x0 <- bp[wi] - bar_width / 2
      x1 <- bp[wi] + bar_width / 2
      y0 <- base; y1 <- base + n_here
      rect(x0, y0, x1, y1, col = ds_colors[all_ds[di]], border = NA)
      stripe_step <- 0.4
      clip(x0, x1, y0, y1)
      span <- (x1 - x0) + (y1 - y0)
      for (s in seq(-span, span, by = stripe_step))
        lines(c(x0, x0 + span), c(y0 + s, y0 + s + span), col = "white", lwd = 1.0)
      do.call("clip", as.list(par("usr")))
      base <- y1
    }
  }

  boundaries <- list(
    list(week = 12, label = "T1 | T2"),
    list(week = 27, label = "T2 | Early Preterm"),
    list(week = 33, label = "Early Preterm | Late Preterm"),
    list(week = 36, label = "Late Preterm | Term")
  )
  bar_spacing <- if (length(bp) > 1) diff(bp[1:2]) else bp[1]
  for (b in boundaries) {
    bc <- which(as.integer(colnames(exact_mat)) == b$week)
    if (length(bc) == 1 && bc < ncol(exact_mat)) {
      x_border <- bp[bc] + bar_spacing / 2
      abline(v = x_border, lty = 2, col = "grey60", lwd = 1.2)
      text(x_border, ymax * 0.97, b$label, cex = 0.5, col = "grey50",
           srt = 90, adj = c(1, 0.5))
    }
  }

  exact_n <- rowSums(exact_mat)
  est_n <- rowSums(est_mat)
  has <- (exact_n + est_n) > 0
  show_ds <- all_ds[has]

  legend("topleft",
         legend = c(sprintf("%s (n=%d)", show_ds, exact_n[show_ds] + est_n[show_ds]),
                    "solid = exact week", "striped = estimated"),
         fill = c(ds_colors[show_ds], "grey60", "grey60"),
         density = c(rep(NA, length(show_ds)), NA, 20),
         angle = c(rep(0, length(show_ds)), 0, 45),
         border = c(rep("white", length(show_ds)), "grey40", "grey40"),
         cex = 0.55, bty = "n", title = "Dataset", ncol = 2)

  dev.off()
  cat("Saved:", plot_path, "\n")

  # --- XLSX ---
  out <- df[, c("sample_id", "dataset", "trimester", "condition",
                "exact_week", "range_weeks", "estimated_week",
                "estimation_method")]
  out$week_source <- ifelse(out$exact_week > 0, "exact", "estimated")
  out <- out[order(out$dataset, out$estimated_week), ]

  wb <- createWorkbook()
  addWorksheet(wb, "gestational_age")
  writeData(wb, "gestational_age", out)

  hs <- createStyle(textDecoration = "bold", fgFill = "#D9E1F2", border = "Bottom")
  addStyle(wb, "gestational_age", hs, rows = 1, cols = 1:ncol(out))

  est_rows <- which(out$week_source == "estimated") + 1
  if (length(est_rows) > 0) {
    addStyle(wb, "gestational_age", createStyle(fgFill = "#FFF3CD"),
             rows = est_rows, cols = 1:ncol(out), gridExpand = TRUE)
  }

  setColWidths(wb, "gestational_age", cols = 1:ncol(out), widths = "auto")

  xlsx_path <- file.path(output_dir, paste0(file_prefix, ".xlsx"))
  saveWorkbook(wb, xlsx_path, overwrite = TRUE)
  cat("Saved:", xlsx_path, "\n")

  cat(sprintf("\nTotal: %d samples (%d exact, %d estimated)\n", nrow(out), n_exact, n_est))
  print(table(out$dataset, out$week_source))
}

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) {
    cat("Usage: Rscript plot_gestational_age_by_tissue.R <config.yaml>\n")
    quit(status = 1)
  }
  plot_gestational_age_by_tissue(args[1])
}
