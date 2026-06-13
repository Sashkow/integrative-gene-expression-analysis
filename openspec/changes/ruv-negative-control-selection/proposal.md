## Why

The existing RUV control gene selection (`build_ruv_control_genes.R`) starts from the Eisenberg-Levanon housekeeping gene list (3,804 symbols) and filters by empirical DE status across multiple pipeline runs. This top-down approach misses potentially stable genes not in the housekeeping catalog. A complementary bottom-up approach — selecting genes empirically from datasets that span both trimesters (GSE9984/Mikheev and GSE100051/Soncin) based on expression stability — would identify candidates grounded in the actual data, then cross-validate against literature and known DEGs.

## What Changes

- New R script that selects RUV negative control candidate genes empirically from Mikheev (GSE9984) and Soncin (GSE100051) datasets using a multi-filter pipeline:
  1. Expression level filter: retain genes in the 25th–90th percentile of average expression (avoids noise floor at low end and probe saturation at high end)
  2. Low variance filter: select genes with lowest coefficient of variation (CV = SD/mean) within the expression band
  3. DE exclusion: remove genes that are differentially expressed (FDR < 0.05) in either dataset's 1st-vs-2nd trimester contrast
  4. Housekeeping cross-reference: flag which candidates overlap with the Eisenberg-Levanon housekeeping list
  5. Literature DEG exclusion: remove genes known to be differentially expressed in 1st-vs-2nd trimester placenta from qPCR/blot studies
- Output: ranked gene list with annotations (expression level, CV, DE status, housekeeping flag, literature verdict), plus a plain ENTREZ ID list for direct use in the RUV pipeline

## Capabilities

### New Capabilities
- `empirical-control-gene-selection`: Bottom-up selection of RUV negative control genes from expression data using expression-level filtering, variance ranking, DE exclusion, and cross-referencing against housekeeping lists and literature DEGs

### Modified Capabilities
<!-- None — this adds a new script alongside the existing build_ruv_control_genes.R -->

## Impact

- New script in `scripts/integrative_analysis/phase2b_direct_merge/` or `scripts/integrative_analysis/article_validation/`
- Reads existing per-dataset expression TSVs from `data/mapped/`
- Reads existing Eisenberg-Levanon housekeeping list from `data/reference/`
- Reads existing pipeline DE results for DE exclusion
- Output to `output/article_validation/ruv_empirical_control_genes/`
- No changes to existing pipeline code — this is an analysis/validation script
