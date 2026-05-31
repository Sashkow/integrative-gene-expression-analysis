## Why

The article currently reports results from a 7-dataset integration (GSE100051, GSE122214, GSE22490, GSE28551, GSE37901, GSE93520, GSE9984). Since then, GSE22490 was identified as having decidual tissue contamination in its protocol, making it unsuitable for a pure placenta/chorionic villi comparison. A cleaner 6-dataset run (`phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko`) has already been completed with updated preprocessing and reference-batch ComBat, yielding 447 DEGs (vs 380 in the 7ds run) from the same 17,531 testable genes on 117 samples. The article text, tables, figures, and supplementary data need to be updated to reflect the 6ds results.

## What Changes

- Replace all references to "seven datasets" / "7 datasets" with "six datasets" / "6 datasets" throughout `main.tex`
- Remove GSE22490 from dataset listings, Table 1, tissue descriptions, and availability statements
- Update key numerical results: DEG counts (380 → 447, up 369 / down 78), intersection DEGs (182 → 221), sample counts (123 → 117), imputation validation correlation (0.823 → 0.816)
- Update normalization method references from `combat` to `combat_ref` (reference-batch ComBat with GSE100051 as reference)
- Replace figures (PCA, NA staircase, validation plots) with 6ds versions from the new output directory
- Update `data/gestational_age_7ds.xlsx` or replace with 6ds equivalent
- Update supplementary data references to point to 6ds output files
- Add brief mention of GSE22490 exclusion rationale (decidual contamination) in Methods

## Capabilities

### New Capabilities

- `article-text-6ds-update`: Update all prose, numbers, and dataset references in `main.tex` from 7ds to 6ds results
- `article-figures-6ds-update`: Replace article figures and data files with 6ds output versions
- `article-supplements-6ds-update`: Update supplementary materials and data availability references

### Modified Capabilities

## Impact

- `articles/imputation_article/main.tex` — extensive text, table, and figure reference changes
- `articles/imputation_article/figures/` — replaced PCA, staircase, and validation plots
- `articles/imputation_article/data/` — gestational age data file replaced
- `articles/imputation_article/supplement1/`, `supplement2/` — updated data tables
- Source data: `output/yehor_sashko/phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko/`
