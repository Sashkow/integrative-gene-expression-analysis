#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(openxlsx))

setwd("/home/shivers/a/r/igea-r/scripts/integrative-gene-expression-analysis")

pheno <- read.csv("data/phenodata/samples_cvs_ga_matched.csv",
                   stringsAsFactors = FALSE, quote = "\"")

cvs_weeks <- c("10", "11", "12")

# ── CVS reference samples ────────────────────────────────────────────────────
cvs <- pheno[pheno$secondaryaccession == "GSE12767" &
             pheno$Diagnosis == "Healthy", ]
cvs_df <- data.frame(
  Sample = cvs$arraydatafile_exprscolumnnames,
  Dataset = "GSE12767",
  GA_week = cvs$Gestational.Age,
  Diagnosis = cvs$Diagnosis,
  Specimen = cvs$Biological.Specimen,
  stringsAsFactors = FALSE
)
cvs_df <- cvs_df[order(as.numeric(cvs_df$GA_week)), ]

# ── GSE100051 ─────────────────────────────────────────────────────────────────
g100 <- pheno[pheno$secondaryaccession == "GSE100051" &
              pheno$Diagnosis == "Healthy" &
              pheno$Gestational.Age %in% cvs_weeks, ]
g100_df <- data.frame(
  Sample = g100$arraydatafile_exprscolumnnames,
  Dataset = "GSE100051",
  Platform = "Illumina HumanHT-12 v4",
  GA_week = g100$Gestational.Age,
  Specimen = g100$Biological.Specimen,
  CVS_at_same_week = sapply(g100$Gestational.Age, function(w)
    sum(cvs_df$GA_week == w)),
  stringsAsFactors = FALSE
)
g100_df <- g100_df[order(as.numeric(g100_df$GA_week)), ]

# ── GSE93520 ──────────────────────────────────────────────────────────────────
g93 <- pheno[pheno$secondaryaccession == "GSE93520" &
             pheno$Diagnosis == "Healthy" &
             pheno$Gestational.Age %in% cvs_weeks, ]
g93_df <- data.frame(
  Sample = g93$arraydatafile_exprscolumnnames,
  Dataset = "GSE93520",
  Platform = "Agilent 4x44K",
  GA_week = g93$Gestational.Age,
  Specimen = g93$Biological.Specimen,
  CVS_at_same_week = sapply(g93$Gestational.Age, function(w)
    sum(cvs_df$GA_week == w)),
  stringsAsFactors = FALSE
)
g93_df <- g93_df[order(as.numeric(g93_df$GA_week)), ]

# ── GSE28551 ──────────────────────────────────────────────────────────────────
g28 <- pheno[pheno$secondaryaccession == "GSE28551" &
             pheno$Diagnosis == "Healthy" &
             pheno$Gestational.Age.Category == "First Trimester", ]
g28_df <- data.frame(
  Sample = g28$arraydatafile_exprscolumnnames,
  Dataset = "GSE28551",
  Platform = "ABI Human Genome v2",
  GA_week = "9-12 (per-sample unknown)",
  Specimen = g28$Biological.Specimen,
  Note = "Group mean 71.2 +/- 8 days (~10.2 weeks)",
  stringsAsFactors = FALSE
)

# ── Not viable datasets ───────────────────────────────────────────────────────
not_viable <- data.frame(
  Dataset = c("GSE22490", "GSE37653", "GSE122214", "GSE9984",
              "GSE55439", "GSE107824"),
  Platform = c("Affymetrix HG-U133Plus2", "NimbleGen", "Affymetrix HG-U133Plus2",
               "Affymetrix HG-U133Plus2", "Illumina HumanHT-12", "Affymetrix HG-U133Plus2"),
  GA_range = c("4-14 (1 at GA 11)", "6-8", "7-8", "6.4-8.4", "unknown", "8-12"),
  n_at_GA_10_12 = c(1, 0, 0, 0, 0, 4),
  Reason = c("Only 1 healthy sample at GA 11",
             "GA 6-8, no overlap with CVS GA 10-12",
             "GA 7-8, no overlap with CVS GA 10-12",
             "GA 6.4-8.4, no overlap with CVS GA 10-12",
             "GA unknown, preservation methods study",
             "Cultured cytotrophoblasts, not tissue; not in mapped data"),
  stringsAsFactors = FALSE
)

# ── Summary sheet ─────────────────────────────────────────────────────────────
summary_df <- data.frame(
  Dataset = c("GSE12767 (CVS)", "GSE100051", "GSE93520", "GSE28551"),
  Platform = c("Affymetrix HG-U133Plus2", "Illumina HumanHT-12 v4",
               "Agilent 4x44K", "ABI Human Genome v2"),
  Role = c("Reference (CVS)", "Termination", "Termination", "Termination"),
  n_first_trim = c(nrow(cvs_df), nrow(g100_df), nrow(g93_df), nrow(g28_df)),
  GA_weeks_matched = c("10, 11, 12", "10, 11, 12", "10", "9-12 (unknown per-sample)"),
  n_at_GA_10 = c(sum(cvs_df$GA_week == "10"), sum(g100_df$GA_week == "10"),
                 sum(g93_df$GA_week == "10"), NA),
  n_at_GA_11 = c(sum(cvs_df$GA_week == "11"), sum(g100_df$GA_week == "11"),
                 sum(g93_df$GA_week == "11"), NA),
  n_at_GA_12 = c(sum(cvs_df$GA_week == "12"), sum(g100_df$GA_week == "12"),
                 sum(g93_df$GA_week == "12"), NA),
  Viable = c("Reference", "Strong", "Usable (GA 10 only)", "Usable (no per-sample GA)"),
  stringsAsFactors = FALSE
)

# ── Write xlsx ────────────────────────────────────────────────────────────────
wb <- createWorkbook()

addWorksheet(wb, "Summary")
writeData(wb, "Summary", summary_df)
freezePane(wb, "Summary", firstRow = TRUE)
setColWidths(wb, "Summary", cols = 1:ncol(summary_df), widths = "auto")

addWorksheet(wb, "CVS_samples")
writeData(wb, "CVS_samples", cvs_df)
freezePane(wb, "CVS_samples", firstRow = TRUE)
setColWidths(wb, "CVS_samples", cols = 1:ncol(cvs_df), widths = "auto")

addWorksheet(wb, "GSE100051_matches")
writeData(wb, "GSE100051_matches", g100_df)
freezePane(wb, "GSE100051_matches", firstRow = TRUE)
setColWidths(wb, "GSE100051_matches", cols = 1:ncol(g100_df), widths = "auto")

addWorksheet(wb, "GSE93520_matches")
writeData(wb, "GSE93520_matches", g93_df)
freezePane(wb, "GSE93520_matches", firstRow = TRUE)
setColWidths(wb, "GSE93520_matches", cols = 1:ncol(g93_df), widths = "auto")

addWorksheet(wb, "GSE28551_first_trim")
writeData(wb, "GSE28551_first_trim", g28_df)
freezePane(wb, "GSE28551_first_trim", firstRow = TRUE)
setColWidths(wb, "GSE28551_first_trim", cols = 1:ncol(g28_df), widths = "auto")

addWorksheet(wb, "Not_viable")
writeData(wb, "Not_viable", not_viable)
freezePane(wb, "Not_viable", firstRow = TRUE)
setColWidths(wb, "Not_viable", cols = 1:ncol(not_viable), widths = "auto")

out_path <- "articles/cvs_vs_termination_compare/ga_week_matches.xlsx"
saveWorkbook(wb, out_path, overwrite = TRUE)
cat(sprintf("Wrote: %s\n", out_path))
cat(sprintf("  CVS: %d samples\n", nrow(cvs_df)))
cat(sprintf("  GSE100051: %d matched\n", nrow(g100_df)))
cat(sprintf("  GSE93520: %d matched\n", nrow(g93_df)))
cat(sprintf("  GSE28551: %d first-trim (no per-sample GA)\n", nrow(g28_df)))
