#!/usr/bin/env Rscript
#' Select negative control genes for RUV (Remove Unwanted Variation).
#'
#' Starting set: Eisenberg & Levanon (2013) human housekeeping genes (3,804
#' symbols, ~3,505 mapped to Entrez).  Gene lists are built separately for
#' the 1st-vs-2nd trimester (1_2) and 2nd-vs-Term (2_3) contrasts.
#'
#' ── Inclusion / exclusion criteria (applied per contrast) ───────────────────
#'
#'   1. EXCLUDE if the gene is a DEG (adj.P.Val < 0.05) in the **balanced**
#'      run (only datasets contributing to both sides of the contrast, ComBat)
#'      OR the **batch-in-limma** run (all datasets, batch as covariate).
#'      Both are conservative approaches: balanced uses only datasets
#'      with samples on both sides of the contrast, and batch-in-limma
#'      avoids ComBat entirely.  If even these flag a gene, it is likely
#'      genuinely DE rather than a batch artefact.
#'
#'   2. EXCLUDE if the gene is **absent from both balanced and batch-in-limma
#'      runs** (filtered out by coverage/variance) AND is a DEG in **every**
#'      other run where it is present (all_datasets, restoration, 2nd_trim_only).
#'      Rationale: no evidence of stability, only evidence of DE.
#'
#'   3. INCLUDE (override empirical DEG status) if placenta-specific qPCR
#'      literature ranks the gene as **stable** — UNLESS our pipeline shows it
#'      is a strong DEG (|logFC| > 0.5 AND FDR < 0.05) in the balanced or
#'      batch-in-limma run.  Literature sources: Meller et al. 2005,
#'      Drewlo et al. 2012, Lanoix et al. 2012, St-Pierre et al. 2017,
#'      Murthi et al. 2008.
#'
#'   4. EXCLUDE (override empirical stability) if literature flags the gene
#'      as **unstable** in placenta (e.g. ACTB, GAPDH, RPLP0).
#'
#'   5. Genes with **no literature verdict** ("neutral" or "none") follow
#'      empirical data only: included if non-DEG, excluded if DEG or absent.
#'
#' ── Outputs ─────────────────────────────────────────────────────────────────
#'   - ruv_control_gene_selection.xlsx   per-gene decisions with reasoning
#'   - ruv_control_genes_1_2.txt         plain Entrez ID list for 1_2
#'   - ruv_control_genes_2_3.txt         plain Entrez ID list for 2_3

suppressPackageStartupMessages({
  library(openxlsx)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

setwd("/home/shivers/a/r/igea-r/scripts/integrative-gene-expression-analysis")

fdr_cutoff       <- 0.05
strong_logfc     <- 0.5

# ── Load Eisenberg-Levanon housekeeping gene list ────────────────────────────
el_file <- "data/reference/integration_methods_references/ruv_housekeeping_genes_search/eisenberg_levanon_HK_genes.txt"
el_raw  <- read.delim(el_file, header = FALSE, stringsAsFactors = FALSE, strip.white = TRUE)
el_symbols <- unique(trimws(el_raw$V1))
cat(sprintf("Eisenberg-Levanon list: %d symbols\n", length(el_symbols)))

el_map <- suppressMessages(
  AnnotationDbi::select(org.Hs.eg.db, keys = el_symbols, columns = "ENTREZID", keytype = "SYMBOL")
)
el_map <- el_map[!is.na(el_map$ENTREZID), ]
el_map <- el_map[!duplicated(el_map$ENTREZID), ]
cat(sprintf("Mapped to Entrez: %d\n", nrow(el_map)))

# ── Literature verdicts for genes in our 32-gene panel ───────────────────────
# verdict categories: "stable", "unstable", "mid-rank" (treated as neutral)
# Per-contrast where possible; otherwise overall
lit_verdicts <- list(
  # Meller 2005 + Drewlo 2012 + Lanoix 2012 + St-Pierre 2017 + Murthi 2008
  YWHAZ  = list(overall = "stable",  note = "Top stable in Meller, Drewlo, St-Pierre; stable in all our contrasts"),
  TBP    = list(overall = "stable",  note = "Top stable in Meller; best pair in St-Pierre geNorm; least stable in Lanoix/Murthi"),
  SDHA   = list(overall = "stable",  note = "Top stable in Meller (rank 2); inconsistent in Lanoix/Murthi"),
  TOP1   = list(overall = "stable",  note = "Top stable in Drewlo (rank 1); stable in Lanoix NormFinder"),
  CYC1   = list(overall = "stable",  note = "Top stable in Drewlo (rank 2)"),
  HPRT1  = list(v_1_2 = "stable", v_2_3 = "neutral", note = "Stable 1st-vs-3rd (Drewlo); top in PE (Lanoix); sex-biased (St-Pierre); strongly DE in 2_3 in our data"),
  IPO8   = list(overall = "stable",  note = "Top NormFinder in St-Pierre; mid-rank in Drewlo"),
  RPL13A = list(overall = "neutral", note = "Sex-biased in St-Pierre; mid-rank in Drewlo"),
  POLR2A = list(overall = "neutral", note = "Mid-rank in Meller and Drewlo"),
  PPIA   = list(overall = "stable",  note = "Top stable in Lanoix (PE+GDM); stable in St-Pierre"),
  ACTB   = list(overall = "unstable", note = "Least stable in Meller and Drewlo; least stable in Murthi"),
  GAPDH  = list(overall = "unstable", note = "Least stable in Meller and Drewlo; context-dependent in Lanoix/Murthi"),
  RPLP0  = list(overall = "unstable", note = "Least stable in Drewlo"),
  B2M    = list(overall = "neutral", note = "Mid-rank in Meller and Drewlo"),
  HMBS   = list(overall = "neutral", note = "Mid-rank in Meller and Drewlo"),
  UBC    = list(overall = "neutral", note = "Mid-rank in Meller, Drewlo, St-Pierre"),
  PGK1   = list(overall = "neutral", note = "Mid-rank in Drewlo and St-Pierre"),
  GUSB   = list(overall = "neutral", note = "Sex-biased in St-Pierre; mid-rank in Drewlo"),
  EIF4A2 = list(overall = "neutral", note = "Mid-rank in Drewlo"),
  NONO   = list(overall = "neutral", note = "Stable in males (St-Pierre)"),
  PSMC4  = list(overall = "neutral", note = "Only in St-Pierre, mid-rank"),
  PUM1   = list(overall = "neutral", note = "Only in St-Pierre, stable"),
  HBB    = list(overall = "neutral", note = "Only in St-Pierre"),
  ALAS1  = list(overall = "neutral", note = "Only in St-Pierre"),
  CDKN1A = list(overall = "neutral", note = "Only in St-Pierre"),
  G6PD   = list(overall = "neutral", note = "Only in St-Pierre"),
  HSP90AB1 = list(overall = "neutral", note = "Only in St-Pierre"),
  LDHA   = list(overall = "neutral", note = "Only in St-Pierre"),
  PPIH   = list(overall = "neutral", note = "Only in St-Pierre"),
  RPL30  = list(overall = "neutral", note = "Best DeltaCq*M pair in St-Pierre"),
  RPS18  = list(overall = "neutral", note = "Sex-biased in St-Pierre"),
  TFRC   = list(overall = "neutral", note = "Only in St-Pierre")
)

get_lit_verdict <- function(symbol, contrast) {
  v <- lit_verdicts[[symbol]]
  if (is.null(v)) return(list(verdict = "none", note = "Not in literature panel"))
  key <- paste0("v_", contrast)
  verdict <- if (!is.null(v[[key]])) v[[key]] else if (!is.null(v$overall)) v$overall else "neutral"
  list(verdict = verdict, note = v$note)
}

# ── Load DE results ──────────────────────────────────────────────────────────
de_files <- list(
  "1_2" = list(
    all_datasets   = "output/phase2b_combat/phase2b_1_2_all_datasets/difexp_none_combat.tsv",
    batch_in_limma = "output/phase2b_batch_in_limma/phase2b_1_2_all_datasets_batch_in_limma/difexp_none_batch_in_limma.tsv",
    balanced       = "output/phase2b_combat/phase2b_1_2_balanced/difexp_none_combat.tsv",
    restoration    = "output/phase2b_combat/phase2b_1_2_restoration_blockmask_imputed_0/difexp_none_combat.tsv",
    second_only    = "output/phase2b_combat/phase2b_1_2_2nd_trim_only/difexp_none_combat.tsv"
  ),
  "2_3" = list(
    all_datasets   = "output/phase2b_combat/phase2b_2_3_all_datasets/difexp_none_combat.tsv",
    batch_in_limma = "output/phase2b_batch_in_limma/phase2b_2_3_all_datasets_batch_in_limma/difexp_none_batch_in_limma.tsv",
    balanced       = "output/phase2b_combat/phase2b_2_3_balanced/difexp_none_combat.tsv",
    restoration    = "output/phase2b_combat/phase2b_2_3_restoration_blockmask_imputed_0/difexp_none_combat.tsv",
    second_only    = "output/phase2b_combat/phase2b_2_3_2nd_trim_only/difexp_none_combat.tsv"
  )
)

load_de <- function(path) {
  if (!file.exists(path)) return(NULL)
  read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
             check.names = FALSE, quote = "")
}

de_tables <- list()
for (contrast in names(de_files)) {
  de_tables[[contrast]] <- list()
  for (run_type in names(de_files[[contrast]])) {
    d <- load_de(de_files[[contrast]][[run_type]])
    if (!is.null(d)) {
      de_tables[[contrast]][[run_type]] <- d
      cat(sprintf("  Loaded %s/%s: %d genes\n", contrast, run_type, nrow(d)))
    }
  }
}

lookup <- function(entrez_id, contrast, run_type) {
  d <- de_tables[[contrast]][[run_type]]
  if (is.null(d)) return(list(logFC = NA, fdr = NA, present = FALSE))
  idx <- which(as.character(d$gene) == entrez_id)
  if (length(idx) == 0) return(list(logFC = NA, fdr = NA, present = FALSE))
  list(logFC = d$logFC[idx[1]], fdr = d$adj.P.Val[idx[1]], present = TRUE)
}

is_deg <- function(res) {
  !is.na(res$fdr) && res$fdr < fdr_cutoff
}

is_strong_deg <- function(res) {
  !is.na(res$fdr) && res$fdr < fdr_cutoff && !is.na(res$logFC) && abs(res$logFC) > strong_logfc
}

# ── Evaluate each gene for each contrast ─────────────────────────────────────
evaluate_gene <- function(entrez_id, symbol, contrast) {
  bal <- lookup(entrez_id, contrast, "balanced")
  bil <- lookup(entrez_id, contrast, "batch_in_limma")
  alld <- lookup(entrez_id, contrast, "all_datasets")
  rest <- lookup(entrez_id, contrast, "restoration")
  sec  <- lookup(entrez_id, contrast, "second_only")

  lit <- get_lit_verdict(symbol, contrast)
  reasons <- character(0)

  # Rule 1: Exclude if DEG in balanced or batch_in_limma
  deg_in_primary <- FALSE
  if (is_deg(bal)) {
    deg_in_primary <- TRUE
    reasons <- c(reasons, sprintf("DEG in balanced (logFC=%.2f, FDR=%.2e)", bal$logFC, bal$fdr))
  }
  if (is_deg(bil)) {
    deg_in_primary <- TRUE
    reasons <- c(reasons, sprintf("DEG in batch_in_limma (logFC=%.2f, FDR=%.2e)", bil$logFC, bil$fdr))
  }

  # Rule 2: If not present in balanced/limma, check if DEG in all runs where present
  not_in_primary <- !bal$present && !bil$present
  if (not_in_primary) {
    other_runs <- list(alld, rest, sec)
    present_runs <- Filter(function(r) r$present, other_runs)
    if (length(present_runs) > 0) {
      deg_in_all_present <- all(sapply(present_runs, is_deg))
      if (deg_in_all_present) {
        deg_in_primary <- TRUE
        reasons <- c(reasons, "Not in balanced/limma; DEG in ALL other runs where present")
      } else {
        reasons <- c(reasons, sprintf("Not in balanced/limma; present in %d other runs, DEG in %d",
                                       length(present_runs),
                                       sum(sapply(present_runs, is_deg))))
      }
    } else {
      reasons <- c(reasons, "Not present in any run")
    }
  }

  # Collect best available data for reporting
  all_results <- list(balanced = bal, batch_in_limma = bil,
                      all_datasets = alld, restoration = rest, second_only = sec)
  present_results <- Filter(function(r) r$present, all_results)
  best_logfc <- NA_real_
  best_fdr   <- NA_real_
  if (length(present_results) > 0) {
    fdrs <- sapply(present_results, function(r) r$fdr)
    best_idx <- which.min(fdrs)
    best_logfc <- present_results[[best_idx]]$logFC
    best_fdr   <- present_results[[best_idx]]$fdr
  }

  # Rule 3: Literature stable override (unless strong DEG)
  if (lit$verdict == "stable") {
    if (deg_in_primary && is_strong_deg(bal) || is_strong_deg(bil)) {
      reasons <- c(reasons,
                   sprintf("Literature: stable (%s) but OVERRIDDEN by strong DEG |logFC|>0.5", lit$note))
      decision <- "EXCLUDE"
    } else if (deg_in_primary) {
      reasons <- c(reasons,
                   sprintf("Literature: stable (%s); weak DEG overridden by literature", lit$note))
      decision <- "INCLUDE"
    } else {
      reasons <- c(reasons, sprintf("Literature: stable (%s)", lit$note))
      decision <- "INCLUDE"
    }
  } else if (lit$verdict == "unstable") {
    # Rule 4: Literature unstable override
    reasons <- c(reasons, sprintf("Literature: UNSTABLE (%s)", lit$note))
    decision <- "EXCLUDE"
  } else {
    # No literature override — use empirical data
    if (deg_in_primary) {
      decision <- "EXCLUDE"
    } else if (not_in_primary && length(present_results) == 0) {
      decision <- "EXCLUDE"
      reasons <- c(reasons, "No empirical data available")
    } else {
      decision <- "INCLUDE"
      if (length(reasons) == 0) reasons <- "Stable in all tested runs"
    }
  }

  list(
    decision   = decision,
    reason     = paste(reasons, collapse = "; "),
    lit_verdict = lit$verdict,
    lit_note    = lit$note,
    logFC_bal   = bal$logFC,
    fdr_bal     = bal$fdr,
    logFC_bil   = bil$logFC,
    fdr_bil     = bil$fdr,
    best_logfc  = best_logfc,
    best_fdr    = best_fdr,
    n_runs_present = sum(sapply(all_results, function(r) r$present)),
    n_runs_deg     = sum(sapply(all_results, is_deg))
  )
}

# ── Build output for both contrasts ──────────────────────────────────────────
build_contrast_table <- function(contrast) {
  rows <- list()
  for (i in seq_len(nrow(el_map))) {
    sym <- el_map$SYMBOL[i]
    eid <- el_map$ENTREZID[i]
    ev  <- evaluate_gene(eid, sym, contrast)
    rows[[i]] <- data.frame(
      symbol          = sym,
      entrez_id       = eid,
      decision        = ev$decision,
      reason          = ev$reason,
      lit_verdict     = ev$lit_verdict,
      lit_note        = ev$lit_note,
      logFC_balanced  = ev$logFC_bal,
      FDR_balanced    = ev$fdr_bal,
      logFC_bil       = ev$logFC_bil,
      FDR_bil         = ev$fdr_bil,
      best_logFC      = ev$best_logfc,
      best_FDR        = ev$best_fdr,
      n_runs_present  = ev$n_runs_present,
      n_runs_deg      = ev$n_runs_deg,
      stringsAsFactors = FALSE
    )
  }
  df <- do.call(rbind, rows)
  df <- df[order(df$decision, df$symbol), ]
  rownames(df) <- NULL
  df
}

cat("\n=== Evaluating 1_2 contrast ===\n")
df_1_2 <- build_contrast_table("1_2")
cat(sprintf("  INCLUDE: %d, EXCLUDE: %d\n",
            sum(df_1_2$decision == "INCLUDE"), sum(df_1_2$decision == "EXCLUDE")))

cat("\n=== Evaluating 2_3 contrast ===\n")
df_2_3 <- build_contrast_table("2_3")
cat(sprintf("  INCLUDE: %d, EXCLUDE: %d\n",
            sum(df_2_3$decision == "INCLUDE"), sum(df_2_3$decision == "EXCLUDE")))

# ── Summary sheet ────────────────────────────────────────────────────────────
overlap_include <- intersect(
  df_1_2$entrez_id[df_1_2$decision == "INCLUDE"],
  df_2_3$entrez_id[df_2_3$decision == "INCLUDE"]
)
cat(sprintf("\nIncluded in BOTH contrasts: %d genes\n", length(overlap_include)))

summary_df <- data.frame(
  metric = c(
    "EL genes mapped to Entrez",
    "1_2: included", "1_2: excluded",
    "2_3: included", "2_3: excluded",
    "Included in BOTH contrasts",
    "Literature overrides (stable -> include despite weak DEG)",
    "Literature overrides (unstable -> exclude despite stable data)"
  ),
  value = c(
    nrow(el_map),
    sum(df_1_2$decision == "INCLUDE"), sum(df_1_2$decision == "EXCLUDE"),
    sum(df_2_3$decision == "INCLUDE"), sum(df_2_3$decision == "EXCLUDE"),
    length(overlap_include),
    sum(df_1_2$decision == "INCLUDE" & df_1_2$lit_verdict == "stable" & df_1_2$n_runs_deg > 0) +
      sum(df_2_3$decision == "INCLUDE" & df_2_3$lit_verdict == "stable" & df_2_3$n_runs_deg > 0),
    sum(df_1_2$decision == "EXCLUDE" & df_1_2$lit_verdict == "unstable") +
      sum(df_2_3$decision == "EXCLUDE" & df_2_3$lit_verdict == "unstable")
  ),
  stringsAsFactors = FALSE
)

# ── Legend ────────────────────────────────────────────────────────────────────
legend_df <- data.frame(
  Column = c(
    "symbol", "entrez_id", "decision", "reason",
    "lit_verdict", "lit_note",
    "logFC_balanced", "FDR_balanced",
    "logFC_bil", "FDR_bil",
    "best_logFC", "best_FDR",
    "n_runs_present", "n_runs_deg"
  ),
  Meaning = c(
    "HGNC gene symbol.",
    "NCBI Entrez Gene ID.",
    "INCLUDE or EXCLUDE for RUV negative control set.",
    "Chain of reasoning: which rules applied and why.",
    "Literature verdict: stable / unstable / neutral / none.",
    "Literature details (study references and findings).",
    "logFC from balanced run (datasets with both contrast groups, ComBat).",
    "adj.P.Val from balanced run.",
    "logFC from batch-in-limma run (all datasets, no ComBat).",
    "adj.P.Val from batch-in-limma run.",
    "logFC from run with lowest FDR (most significant across all runs).",
    "Lowest FDR across all runs (worst case).",
    "Number of pipeline runs (out of 5) where this gene passed filters.",
    "Number of runs where gene is DEG (FDR < 0.05)."
  ),
  stringsAsFactors = FALSE
)

# ── Write XLSX ───────────────────────────────────────────────────────────────
wb <- createWorkbook()

red_style   <- createStyle(fontColour = "#CC0000", bgFill = "#FFE0E0")
green_style <- createStyle(fontColour = "#006600", bgFill = "#E0FFE0")
incl_style  <- createStyle(bgFill = "#E0FFE0")
excl_style  <- createStyle(bgFill = "#FFE0E0")

write_contrast_sheet <- function(wb, sheet_name, df) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, df, withFilter = TRUE)
  freezePane(wb, sheet_name, firstActiveRow = 2, firstActiveCol = 4)
  setColWidths(wb, sheet_name, cols = seq_len(ncol(df)), widths = "auto")
  setColWidths(wb, sheet_name, cols = which(colnames(df) == "reason"), widths = 80)
  setColWidths(wb, sheet_name, cols = which(colnames(df) == "lit_note"), widths = 50)

  dec_col <- which(colnames(df) == "decision")
  conditionalFormatting(wb, sheet_name, cols = dec_col, rows = 2:(nrow(df) + 1),
                        rule = '=="INCLUDE"', style = incl_style)
  conditionalFormatting(wb, sheet_name, cols = dec_col, rows = 2:(nrow(df) + 1),
                        rule = '=="EXCLUDE"', style = excl_style)

  for (col_name in c("FDR_balanced", "FDR_bil", "best_FDR")) {
    col_idx <- which(colnames(df) == col_name)
    if (length(col_idx) == 1) {
      conditionalFormatting(wb, sheet_name, cols = col_idx, rows = 2:(nrow(df) + 1),
                            rule = "<0.05", style = red_style)
    }
  }
}

write_contrast_sheet(wb, "RUV_controls_1_2", df_1_2)
write_contrast_sheet(wb, "RUV_controls_2_3", df_2_3)

addWorksheet(wb, "summary")
writeData(wb, "summary", summary_df)
setColWidths(wb, "summary", cols = 1:2, widths = c(55, 15))

addWorksheet(wb, "legend")
writeData(wb, "legend", legend_df)
setColWidths(wb, "legend", cols = 1:2, widths = c(20, 80))

# ── Also write plain gene lists for direct use ───────────────────────────────
genes_1_2 <- df_1_2$entrez_id[df_1_2$decision == "INCLUDE"]
genes_2_3 <- df_2_3$entrez_id[df_2_3$decision == "INCLUDE"]

writeLines(genes_1_2, "articles/imputation_article/ruv_control_genes_1_2.txt")
writeLines(genes_2_3, "articles/imputation_article/ruv_control_genes_2_3.txt")

out_path <- "articles/imputation_article/ruv_control_gene_selection.xlsx"
saveWorkbook(wb, out_path, overwrite = TRUE)

cat(sprintf("\nWrote: %s\n", out_path))
cat(sprintf("Wrote: ruv_control_genes_1_2.txt (%d genes)\n", length(genes_1_2)))
cat(sprintf("Wrote: ruv_control_genes_2_3.txt (%d genes)\n", length(genes_2_3)))
