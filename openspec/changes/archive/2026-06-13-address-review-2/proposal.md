## Why

Reviewer 2 returned a "minor revision" decision with ~20 findings across the main article, figures, and two supplements. All major concerns from Review 1 were resolved; the remaining items are one must-fix citation error, several minor text/methods improvements, and presentation suggestions. Addressing these completes the revision cycle.

## What Changes

**Main text (items 2, 4, 5, 10):**
- Investigate and describe the 2.41 max-deviation outlier gene in ComBat sensitivity analysis
- Add block-mask RMSE calibration as percentage of expression SD (~53%) alongside existing random-cell calibration
- Break up the long Conclusions paragraph (lines 367-376) into two sentences separating primary results from extended validation
- Add Methods paragraph describing worst-case block mask and hard block mask strategies (currently only in Results)

**Figures (items 13, 15):**
- Increase Figure 2 subplot titles (before/after ComBat 2x2)
- Round r=0.999134 to r=0.999 in Figure 4 scatter title

**Presentation (items 3/17, 6, 7, 8):**
- Move Table 4 (extended block-mask) to supplementary material
- Add brief note on the 20 gained genes absent from Affymetrix platforms in concordance section
- Shorten deconvolution limitation bullet
- Standardise to American English spelling (colour → color)

**Footnote verification (item 9):**
- Verify all `\srcm`/`\srct` footnotes resolve on the branch referenced by `\ghb`

**Supplement 1 (items S1.1-S1.4):**
- Note that main article adopted recommendation #5 (ref-batch ComBat) and reference r=0.9991 result
- Acknowledge HarmonizR as future work (not benchmarked)
- Justify M1 vs M2 imputation choice or acknowledge as limitation
- Remove unused GEO dataset citations

**Supplement 2 (items S2.1-S2.4):**
- **Must-fix**: Correct Lanoix citation author duplication (ref 15)
- Update "14 datasets" / "227 DEGs" references to match current 6ds configuration
- Note which RUV improvement strategies are planned for follow-up
- Identify the 10 RUVinv DEGs

## Capabilities

### New Capabilities
- `main-text-revisions`: Text corrections and additions to the main article (Limitations, Methods, Conclusions, concordance section, ComBat sensitivity)
- `figure-adjustments`: Visual fixes to Figures 2 and 4
- `table-restructure`: Move Table 4 to supplement, adjust table numbering
- `supplement-corrections`: Fixes to Supplements 1 and 2 (citation, consistency, justifications)
- `footnote-verification`: Verify all source-code footnotes resolve on the article branch

### Modified Capabilities

## Impact

- `articles/imputation_article/main.tex` — main text edits, table restructuring, figure adjustments
- `articles/imputation_article/supplements/` — supplement 1 and 2 edits
- R scripts generating Figures 2 and 4 — subplot title size, correlation rounding
- No code logic changes; no pipeline reruns needed
- No new dependencies
