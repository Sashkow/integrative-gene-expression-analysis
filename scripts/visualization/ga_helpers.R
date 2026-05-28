GA_TRI_MAP <- c(
  "First Trimester" = "First trimester",
  "Second Trimester" = "Second trimester",
  "Third Trimester" = "Third trimester",
  "Term" = "Term",
  "Late Preterm" = "Late preterm",
  "Early Preterm" = "Early preterm"
)

GA_RANGES <- list(
  "First trimester"  = "4,5,6,7,8,9,10,11,12",
  "Second trimester" = "13,14,15,16,17,18,19,20,21,22,23,24,25,26,27",
  "Third trimester"  = "37,38,39,40,41",
  "Term"             = "37,38,39,40,41",
  "Late preterm"     = "34,35,36",
  "Early preterm"    = "28,29,30,31,32,33"
)

GA_BASE_PALETTE <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#999999",
  "#882255", "#44AA99", "#332288", "#DDCC77",
  "#117733", "#88CCEE", "#AA4499", "#661100",
  "#6699CC", "#000000", "#E6AB02", "#1B9E77",
  "#AE76A3", "#A65628", "#F781BF", "#66C2A5",
  "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
  "#FFD92F", "#B3B3B3"
)

middle_out <- function(opts) {
  n <- length(opts)
  if (n <= 1) return(opts)
  mid <- ceiling(n / 2)
  result <- integer(n)
  left <- mid; right <- mid + 1; pos <- 1
  result[pos] <- opts[mid]
  while (pos < n) {
    if (right <= n) { pos <- pos + 1; result[pos] <- opts[right]; right <- right + 1 }
    if (pos < n && left > 1) { left <- left - 1; pos <- pos + 1; result[pos] <- opts[left] }
  }
  result
}

parse_ga_rows <- function(pheno, extra_cols = character(0)) {
  rows <- list()
  for (i in seq_len(nrow(pheno))) {
    sid <- pheno$arraydatafile_exprscolumnnames[i]
    ds <- pheno$secondaryaccession[i]
    tri <- unname(GA_TRI_MAP[pheno$Gestational.Age.Category[i]])
    if (is.na(tri)) tri <- pheno$Gestational.Age.Category[i]

    ga_raw <- as.character(pheno$Gestational.Age[i])
    ga_val <- suppressWarnings(as.integer(ga_raw))
    if (is.na(ga_val) || ga_val <= 0) {
      wk <- regmatches(ga_raw, regexpr("\\d+", ga_raw))
      ga_val <- if (length(wk) == 1) as.integer(wk) else 0L
    }

    ga_lower <- suppressWarnings(as.integer(pheno$Gestational.Age.Lower.Bound[i]))
    ga_upper <- suppressWarnings(as.integer(pheno$Gestational.Age.Upper.Bound[i]))
    rw <- ""
    if (ga_val == 0L && !is.na(ga_lower) && !is.na(ga_upper) && ga_lower > 0 && ga_upper > 0)
      rw <- paste(seq(ga_lower, ga_upper), collapse = ",")

    rep_block <- as.character(pheno$technical_replicate_block[i])
    if (is.na(rep_block) || rep_block %in% c("", "_")) rep_block <- NA_character_

    row <- data.frame(
      sample_id = sid, dataset = ds, trimester = tri,
      exact_week = ga_val, range_weeks = rw, rep_block = rep_block,
      stringsAsFactors = FALSE)

    for (nm in names(extra_cols))
      row[[nm]] <- pheno[[extra_cols[[nm]]]][i]

    rows[[length(rows) + 1]] <- row
  }
  do.call(rbind, rows)
}

collapse_rep_blocks <- function(df) {
  df <- df[!duplicated(df$sample_id), ]
  has_block <- !is.na(df$rep_block)
  if (any(has_block)) {
    n_before <- nrow(df)
    keep <- !has_block | !duplicated(df$rep_block)
    df <- df[keep, ]
    cat(sprintf("Collapsed technical replicates: %d -> %d unique samples\n",
                n_before, nrow(df)))
  }
  df
}

fill_range_from_category <- function(df) {
  for (i in which(df$exact_week == 0 & (is.na(df$range_weeks) | df$range_weeks == ""))) {
    tri <- df$trimester[i]
    if (tri %in% names(GA_RANGES))
      df$range_weeks[i] <- GA_RANGES[[tri]]
  }
  df
}

estimate_weeks <- function(df) {
  df$estimated_week <- df$exact_week
  df$estimation_method <- ifelse(df$exact_week > 0, "exact", NA)
  df$range_key <- NA_character_

  range_idx <- which(df$exact_week == 0)
  for (i in range_idx) {
    rw <- df$range_weeks[i]
    if (!is.na(rw) && rw != "" && rw != "0") {
      opts <- as.integer(unlist(strsplit(as.character(rw), ",")))
      opts <- opts[!is.na(opts)]
      df$estimation_method[i] <- "even spread across range (middle-out)"
    } else {
      tri <- df$trimester[i]
      if (tri %in% names(GA_RANGES)) {
        opts <- as.integer(unlist(strsplit(GA_RANGES[[tri]], ",")))
        df$estimation_method[i] <- "even spread across trimester (middle-out)"
      } else {
        next
      }
    }
    df$range_key[i] <- paste(sort(opts), collapse = ",")
  }

  for (key in unique(paste(df$dataset[range_idx], df$range_key[range_idx], sep = "|"))) {
    if (is.na(key) || key == "NA|NA") next
    parts <- strsplit(key, "\\|")[[1]]
    ds_id <- parts[1]; rk <- parts[2]
    if (is.na(rk)) next
    week_options <- as.integer(unlist(strsplit(rk, ",")))
    mo <- middle_out(week_options)
    idx <- which(df$dataset == ds_id & !is.na(df$range_key) &
                 df$range_key == rk & df$exact_week == 0)
    if (length(idx) > 0)
      df$estimated_week[idx] <- rep_len(mo, length(idx))
  }

  df[df$estimated_week > 0, ]
}

make_count_matrix <- function(group_col, week_col, all_groups, all_weeks) {
  wk_labels <- as.character(all_weeks)
  mat <- matrix(0, nrow = length(all_groups), ncol = length(all_weeks),
                dimnames = list(all_groups, wk_labels))
  for (i in seq_along(group_col))
    mat[group_col[i], as.character(week_col[i])] <-
      mat[group_col[i], as.character(week_col[i])] + 1
  mat
}
