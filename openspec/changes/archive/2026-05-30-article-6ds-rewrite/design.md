## Context

The article (`articles/imputation_article/main.tex`) currently describes a 7-dataset integration using output from `output/phase2b_combat/phase2b_1_2_all_datasets/`. A cleaner 6-dataset run already exists at `output/yehor_sashko/phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko/` with reference-batch ComBat (`combat_ref`, ref=GSE100051) and sex covariate in limma. The 6ds run excludes GSE22490 (decidual contamination) and uses Yehor's updated preprocessing for GSE28551.

Key number changes (7ds → 6ds):
- Datasets: 7 → 6 (GSE22490 removed)
- Samples: 123 → 117
- Testable genes (softImpute): 17,531 (unchanged, coincidentally same)
- Intersection genes: 8,260 (unchanged — GSE22490 was same platform as 3 others)
- DEGs (softImpute+combat_ref): 380 → 447 (up 369, down 78)
- DEGs (intersection, combat_ref): 182 → 221 (up 174, down 47)
- Imputation validation r: 0.823 → 0.816
- Normalization label: combat → combat_ref
- Vendors/platforms: 4/4 (unchanged — GSE22490 shared Affy platform)

## Goals / Non-Goals

**Goals:**
- Update all numerical results in main.tex to match 6ds combat_ref output
- Replace all figures with 6ds versions
- Update data files (gestational age xlsx)
- Update supplementary materials
- Add GSE22490 exclusion rationale in Methods
- Ensure internal consistency (no stale 7ds numbers remain)

**Non-Goals:**
- Re-running the pipeline (6ds output already exists)
- Changing the article structure or adding new sections
- Updating the `article-biological-validation` change artifacts (those will need separate update later)
- Rebuilding supplement PDFs (just update .tex source; build is a separate step)

## Decisions

1. **Use `softimpute_combat_ref` as the primary method** — The 6ds config uses reference-batch ComBat (ref=GSE100051). This is the `softimpute_combat_ref` row in method_comparison.csv (447 DEGs). The article should present this as the main result, with `none_combat_ref` (221 DEGs) as the intersection-only baseline.

2. **Source directory** — All data comes from `output/yehor_sashko/phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko/`. Figures are copied from there into `articles/imputation_article/figures/`.

3. **Text update approach** — Systematic find-and-replace for "seven" → "six", "7 datasets" → "6 datasets", then manual review of every numerical value (DEG counts, sample counts, correlation values) against the 6ds summary.txt and method_comparison.csv.

4. **GSE22490 exclusion** — Add one sentence in the Methods dataset-selection paragraph explaining the exclusion due to decidual tissue contamination in the protocol. Remove it from Table 1 and all dataset listings.

## Risks / Trade-offs

- [Stale numbers] Some numerical values may be embedded in prose that's hard to grep for → Mitigation: systematic line-by-line review of Results and Methods sections after bulk replacements
- [Validation plots] The 6ds output has `combat_ref` variant plots; need to verify figure filenames match LaTeX references → Mitigation: check each `\includegraphics` path after copying
- [Cross-reference with biological-validation change] The open `article-biological-validation` change uses 7ds DEG counts (380, 182) → Mitigation: marked as non-goal; that change's tasks will need updating separately
