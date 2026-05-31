#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(openxlsx)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

setwd("/home/shivers/a/r/igea-r/scripts/integrative-gene-expression-analysis")

# ── Literature data ──────────────────────────────────────────────────────────
# Each study's genes, context, and per-gene verdicts
studies <- list(
  list(
    id          = "Meller2005",
    citation    = "Meller et al. 2005, Placenta 26:S72-S77",
    context     = "All 3 trimesters (1st, 2nd, 3rd); qPCR + geNorm",
    n_samples   = "28 placentas across 3 trimesters",
    genes_tested = c("ACTB","B2M","GAPDH","HMBS","HPRT1","POLR2A","SDHA","TBP","UBC","YWHAZ"),
    verdicts    = list(
      TBP    = "Top stable (rank 1)",
      SDHA   = "Top stable (rank 2)",
      YWHAZ  = "Top stable (rank 3)",
      HPRT1  = "Mid-rank",
      POLR2A = "Mid-rank",
      HMBS   = "Mid-rank",
      UBC    = "Mid-rank",
      B2M    = "Mid-rank",
      ACTB   = "Least stable",
      GAPDH  = "Least stable"
    )
  ),
  list(
    id          = "Drewlo2012",
    citation    = "Drewlo et al. 2012, Placenta 33:893-898",
    context     = "9 clinical groups (gestational ages + pathologies); geNorm + NormFinder",
    n_samples   = "Multiple clinical groups including trimesters",
    genes_tested = c("ACTB","B2M","CYC1","EIF4A2","GAPDH","GUSB","HMBS","HPRT1",
                     "IPO8","PGK1","PPIA","POLR2A","RPL13A","RPLP0","SDHA","TBP",
                     "TOP1","UBC","YWHAZ"),
    verdicts    = list(
      TOP1   = "Top stable (rank 1)",
      CYC1   = "Top stable (rank 2)",
      YWHAZ  = "Top stable (rank 3)",
      HPRT1  = "Stable (validated 1st vs 3rd trim)",
      IPO8   = "Mid-rank",
      SDHA   = "Mid-rank",
      TBP    = "Mid-rank",
      RPL13A = "Mid-rank",
      EIF4A2 = "Mid-rank",
      UBC    = "Mid-rank",
      B2M    = "Mid-rank",
      PGK1   = "Mid-rank",
      PPIA   = "Mid-rank",
      HMBS   = "Mid-rank",
      POLR2A = "Mid-rank",
      GUSB   = "Mid-rank",
      ACTB   = "Least stable",
      GAPDH  = "Least stable",
      RPLP0  = "Least stable"
    )
  ),
  list(
    id          = "Lanoix2012",
    citation    = "Lanoix et al. 2012, Mol Biotechnol 53:61-70",
    context     = "Normal vs PE and Normal vs GDM (term placentas); geNorm + NormFinder",
    n_samples   = "11 normotensive + 11 PE; 8 normal + 8 GDM",
    genes_tested = c("ACTB","GAPDH","HPRT1","PPIA","SDHA","TBP","TOP1","YWHAZ"),
    verdicts    = list(
      HPRT1  = "Top stable (PE: rank 1 geNorm+NormFinder)",
      PPIA   = "Top stable (PE: rank 2; GDM: rank 1-2)",
      TOP1   = "Stable (PE: rank 4 NormFinder)",
      ACTB   = "Stable (PE context only)",
      YWHAZ  = "Inconsistent between tools (PE); stable in GDM",
      SDHA   = "Inconsistent (GDM: rank 1 NormFinder, low geNorm)",
      GAPDH  = "Mid-rank (GDM: top 3 geNorm; PE: unstable)",
      TBP    = "Least stable (both PE and GDM)"
    )
  ),
  list(
    id          = "StPierre2017",
    citation    = "St-Pierre et al. 2017, Sci Rep 7:16923",
    context     = "Male vs female term placentas; 28 genes; geNorm + NormFinder",
    n_samples   = "20 term placentas (10 male, 10 female)",
    genes_tested = c("ACTB","ALAS1","B2M","CDKN1A","G6PD","GAPDH","GUSB","HBB",
                     "HMBS","HPRT1","HSP90AB1","IPO8","LDHA","NONO","PGK1","PPIA",
                     "PPIH","PSMC4","PUM1","RPL13A","RPL30","RPLP0","RPS18","SDHA",
                     "TBP","TFRC","UBC","YWHAZ"),
    verdicts    = list(
      TBP     = "Top stable (geNorm rank 1 pair; NormFinder rank 2)",
      YWHAZ   = "Top stable (geNorm rank 1 pair; NormFinder low in PE)",
      IPO8    = "Top stable (NormFinder rank 1)",
      NONO    = "Stable (male placentas rank 1)",
      PPIA    = "Stable",
      PUM1    = "Stable",
      SDHA    = "Stable",
      RPL30   = "Stable (best novel DeltaCq*M pair with GAPDH)",
      GAPDH   = "Stable (best novel DeltaCq*M pair with RPL30)",
      ACTB    = "Mid-rank",
      B2M     = "Mid-rank",
      HMBS    = "Mid-rank",
      POLR2A  = "Mid-rank",  # not in this study actually
      ALAS1   = "Mid-rank",
      CDKN1A  = "Mid-rank",
      G6PD    = "Mid-rank",
      HSP90AB1= "Mid-rank",
      LDHA    = "Mid-rank",
      PGK1    = "Mid-rank",
      PPIH    = "Mid-rank",
      PSMC4   = "Mid-rank",
      RPL13A  = "Sex-biased (lower Cq in males, p=0.018)",
      RPLP0   = "Mid-rank",
      RPS18   = "Sex-biased (lower Cq in males)",
      TFRC    = "Mid-rank",
      UBC     = "Mid-rank",
      GUSB    = "Sex-biased (lower Cq in males, p=0.017)",
      HPRT1   = "Sex-biased (lower Cq in males, p=0.017; X-linked)",
      HBB     = "Mid-rank"
    )
  ),
  list(
    id          = "Murthi2008",
    citation    = "Murthi et al. 2008, Placenta 29:798-801",
    context     = "Normal vs FGR (fetal growth restriction); 6 genes",
    n_samples   = "Normal + FGR placentas",
    genes_tested = c("ACTB","GAPDH","SDHA","TBP","YWHAZ"),
    verdicts    = list(
      GAPDH  = "Top stable (with 18S and YWHAZ)",
      YWHAZ  = "Top stable",
      SDHA   = "Mid-rank",
      TBP    = "Least stable (with ACTB)",
      ACTB   = "Least stable"
    )
  )
)

# ── Collect all unique gene symbols ──────────────────────────────────────────
all_symbols <- unique(unlist(lapply(studies, function(s) s$genes_tested)))
all_symbols <- sort(all_symbols)
cat(sprintf("Total unique genes across all studies: %d\n", length(all_symbols)))

# ── Map to Entrez IDs ────────────────────────────────────────────────────────
sym_map <- AnnotationDbi::select(org.Hs.eg.db,
                                  keys = all_symbols,
                                  columns = "ENTREZID",
                                  keytype = "SYMBOL")
sym_to_entrez <- setNames(sym_map$ENTREZID, sym_map$SYMBOL)

# ── Load DE results from pipeline runs ───────────────────────────────────────
de_files <- list(
  "1_2_all_datasets"             = "output/phase2b_combat/phase2b_1_2_all_datasets/difexp_none_combat.tsv",
  "2_3_all_datasets"             = "output/phase2b_combat/phase2b_2_3_all_datasets/difexp_none_combat.tsv",
  "1_2_batch_in_limma"           = "output/phase2b_batch_in_limma/phase2b_1_2_all_datasets_batch_in_limma/difexp_none_batch_in_limma.tsv",
  "2_3_batch_in_limma"           = "output/phase2b_batch_in_limma/phase2b_2_3_all_datasets_batch_in_limma/difexp_none_batch_in_limma.tsv",
  "1_2_balanced"                 = "output/phase2b_combat/phase2b_1_2_balanced/difexp_none_combat.tsv",
  "2_3_balanced"                 = "output/phase2b_combat/phase2b_2_3_balanced/difexp_none_combat.tsv"
)

de_data <- list()
for (run_id in names(de_files)) {
  fpath <- de_files[[run_id]]
  if (!file.exists(fpath)) {
    cat(sprintf("  MISSING: %s\n", fpath))
    next
  }
  d <- read.table(fpath, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   check.names = FALSE, quote = "")
  de_data[[run_id]] <- d
  cat(sprintf("  Loaded %s: %d genes\n", run_id, nrow(d)))
}

lookup_de <- function(entrez_id, run_id) {
  if (is.na(entrez_id) || !(run_id %in% names(de_data))) return(list(logFC = NA, fdr = NA))
  d <- de_data[[run_id]]
  idx <- which(as.character(d$gene) == entrez_id)
  if (length(idx) == 0) return(list(logFC = NA, fdr = NA))
  list(logFC = d$logFC[idx[1]], fdr = d$adj.P.Val[idx[1]])
}

# ── Build main table ─────────────────────────────────────────────────────────
rows <- list()
for (sym in all_symbols) {
  eid <- sym_to_entrez[sym]

  tested_in <- character(0)
  verdict_parts <- character(0)
  for (s in studies) {
    if (sym %in% s$genes_tested) {
      tested_in <- c(tested_in, s$id)
      v <- s$verdicts[[sym]]
      if (!is.null(v)) {
        verdict_parts <- c(verdict_parts, paste0(s$id, ": ", v))
      }
    }
  }

  de_12_all <- lookup_de(eid, "1_2_all_datasets")
  de_23_all <- lookup_de(eid, "2_3_all_datasets")
  de_12_bil <- lookup_de(eid, "1_2_batch_in_limma")
  de_23_bil <- lookup_de(eid, "2_3_batch_in_limma")
  de_12_bal <- lookup_de(eid, "1_2_balanced")
  de_23_bal <- lookup_de(eid, "2_3_balanced")

  all_de <- list(de_12_all, de_23_all, de_12_bil, de_23_bil, de_12_bal, de_23_bal)
  is_deg_any <- FALSE
  n_tested <- sum(sapply(all_de, function(x) !is.na(x$fdr)))
  for (de_res in all_de) {
    if (!is.na(de_res$fdr) && de_res$fdr < 0.05) { is_deg_any <- TRUE; break }
  }

  consensus <- if (n_tested == 0) {
    "NOT IN PIPELINE (filtered out in all runs)"
  } else if (is_deg_any) {
    sprintf("UNSTABLE (DEG in >=1 contrast; tested in %d/6 runs)", n_tested)
  } else {
    sprintf("Stable across all tested contrasts (%d/6 runs)", n_tested)
  }

  rows[[length(rows) + 1]] <- data.frame(
    symbol            = sym,
    entrez_id         = ifelse(is.na(eid), "", eid),
    n_studies         = length(tested_in),
    tested_in         = paste(tested_in, collapse = "; "),
    literature_verdicts = paste(verdict_parts, collapse = "\n"),

    logFC_1_2_all     = de_12_all$logFC,
    FDR_1_2_all       = de_12_all$fdr,
    logFC_2_3_all     = de_23_all$logFC,
    FDR_2_3_all       = de_23_all$fdr,
    logFC_1_2_bil     = de_12_bil$logFC,
    FDR_1_2_bil       = de_12_bil$fdr,
    logFC_2_3_bil     = de_23_bil$logFC,
    FDR_2_3_bil       = de_23_bil$fdr,
    logFC_1_2_bal     = de_12_bal$logFC,
    FDR_1_2_bal       = de_12_bal$fdr,
    logFC_2_3_bal     = de_23_bal$logFC,
    FDR_2_3_bal       = de_23_bal$fdr,

    empirical_verdict = consensus,
    stringsAsFactors  = FALSE
  )
}

main_df <- do.call(rbind, rows)
rownames(main_df) <- NULL

# ── Build studies sheet ──────────────────────────────────────────────────────
studies_df <- data.frame(
  study_id    = sapply(studies, `[[`, "id"),
  citation    = sapply(studies, `[[`, "citation"),
  context     = sapply(studies, `[[`, "context"),
  n_samples   = sapply(studies, `[[`, "n_samples"),
  n_genes     = sapply(studies, function(s) length(s$genes_tested)),
  genes       = sapply(studies, function(s) paste(s$genes_tested, collapse = ", ")),
  stringsAsFactors = FALSE
)

# ── Build legend sheet ───────────────────────────────────────────────────────
legend_df <- data.frame(
  Column = c(
    "symbol", "entrez_id", "n_studies", "tested_in", "literature_verdicts",
    "logFC_1_2_all / FDR_1_2_all",
    "logFC_2_3_all / FDR_2_3_all",
    "logFC_1_2_bil / FDR_1_2_bil",
    "logFC_2_3_bil / FDR_2_3_bil",
    "logFC_1_2_bal / FDR_1_2_bal",
    "logFC_2_3_bal / FDR_2_3_bal",
    "empirical_verdict"
  ),
  Meaning = c(
    "HGNC gene symbol.",
    "NCBI Entrez Gene ID.",
    "Number of literature studies testing this gene.",
    "Study IDs that tested this gene.",
    "Per-study stability verdict (rank and notes).",
    "1st vs 2nd Trimester: all 14 datasets, ComBat normalization, no imputation.",
    "2nd Trimester vs Term: all 14 datasets, ComBat normalization, no imputation.",
    "1st vs 2nd Trimester: all 14 datasets, batch-in-limma (no ComBat), no imputation.",
    "2nd Trimester vs Term: all 14 datasets, batch-in-limma (no ComBat), no imputation.",
    "1st vs 2nd Trimester: balanced datasets only (datasets with both groups), ComBat.",
    "2nd Trimester vs Term: balanced datasets only (datasets with both groups), ComBat.",
    "Overall: UNSTABLE if DEG (FDR<0.05) in any contrast/run; Stable otherwise."
  ),
  stringsAsFactors = FALSE
)

# ── Write XLSX ───────────────────────────────────────────────────────────────
wb <- createWorkbook()

# Main sheet
addWorksheet(wb, "reference_genes")
writeData(wb, "reference_genes", main_df, withFilter = TRUE)
freezePane(wb, "reference_genes", firstActiveRow = 2, firstActiveCol = 6)
setColWidths(wb, "reference_genes", cols = seq_len(ncol(main_df)), widths = "auto")

# Conditional formatting: red FDR < 0.05
fdr_cols <- grep("^FDR_", colnames(main_df))
red_style  <- createStyle(fontColour = "#CC0000", bgFill = "#FFE0E0")
green_style <- createStyle(fontColour = "#006600", bgFill = "#E0FFE0")
for (col_idx in fdr_cols) {
  conditionalFormatting(wb, "reference_genes", cols = col_idx, rows = 2:(nrow(main_df) + 1),
                        rule = "<0.05", style = red_style)
  conditionalFormatting(wb, "reference_genes", cols = col_idx, rows = 2:(nrow(main_df) + 1),
                        rule = ">=0.05", style = green_style)
}

# Verdict column styling
verdict_col <- which(colnames(main_df) == "empirical_verdict")
conditionalFormatting(wb, "reference_genes", cols = verdict_col, rows = 2:(nrow(main_df) + 1),
                      rule = 'SEARCH("UNSTABLE", $R2)', style = red_style, type = "expression")

# Studies sheet
addWorksheet(wb, "studies")
writeData(wb, "studies", studies_df)
setColWidths(wb, "studies", cols = 1:6, widths = c(16, 45, 60, 40, 10, 80))

# Legend sheet
addWorksheet(wb, "legend")
writeData(wb, "legend", legend_df)
setColWidths(wb, "legend", cols = 1:2, widths = c(30, 90))

out_path <- "articles/imputation_article/refgene_literature_vs_pipeline.xlsx"
dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
saveWorkbook(wb, out_path, overwrite = TRUE)

cat(sprintf("\nWrote: %s\n", out_path))
cat(sprintf("  %d genes, %d columns\n", nrow(main_df), ncol(main_df)))
cat(sprintf("  Stable: %d, Unstable: %d\n",
            sum(grepl("Stable", main_df$empirical_verdict)),
            sum(grepl("UNSTABLE", main_df$empirical_verdict))))
