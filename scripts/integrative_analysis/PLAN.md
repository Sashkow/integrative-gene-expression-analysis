# Integrative Analysis Improvement Plan

**Last updated**: 2026-05-08

## Problem Statement

Current pipeline loses ~48% of genes through INNER JOIN when merging 8 datasets.
This is unacceptable information loss that likely excludes biologically important genes.

---

## Phase 1: Tiered Gene Coverage Analysis [DONE]

**Goal:** Replace binary "all or nothing" gene intersection with flexible thresholds.

| Tier | Threshold | Use Case |
|------|-----------|----------|
| Tier 1 | 100% (8/8 studies) | High-confidence biomarker discovery |
| Tier 2 | >=75% (6/8 studies) | Balanced primary analysis |
| Tier 3 | >=50% (4/8 studies) | Exploratory analysis |

**Tasks:**
- [x] Add gene coverage reporting to current pipeline
- [x] Implement threshold parameter in merging functions
- [x] Compare DEG results across tiers
- [x] Gene presence matrix exported (`articles/imputation_article/phase2b_gene_presence.xlsx`)

**Results:** Implemented in `phase2b_direct_merge/run_phase2b.R`. Coverage tiers used across all runs. See `phase2b_runs_summary.xlsx` for comprehensive comparison.

---

## Phase 2: Meta-Analysis Methods [DONE]

Implemented in `phase2_meta_analysis/`. Effect-size meta-analysis (DExMA) and rank-based (RankProd) tested. See `research/meta_analysis_vs_merging_comparison.md` for comparison with direct merge.

**Conclusion:** Direct merge with imputation (Phase 2B) was chosen over meta-analysis for the imputation article, as it yielded better gene recovery and allowed fetal sex stratification.

---

## Phase 2B: Direct Merge with Imputation [DONE - Published]

**Published as:** Yehor conference abstract (2026) and imputation article (`articles/imputation_article/`).

### What Was Done

**Imputation methods tested** (via `build_runs_summary_xlsx.R`):
- softImpute (nuclear norm matrix completion) - **selected**
- kNN imputation
- missMDA (PCA-based)
- sample-kNN
- none (baseline)

**Batch correction methods tested:**
- ComBat - **selected**
- DWD - rejected (too aggressive, eliminates biological signal)
- RUV (Remove Unwanted Variation) - tested with k=2,3,4 and variants (bruv, ruvinv)
- Batch as limma covariate (no explicit correction)

**Runs completed** (14 configurations in `phase2b_runs_summary.xlsx`):
- 1_2 and 2_3 cohorts x {all_datasets, balanced, 2nd_trim_only, restoration_blockmask}
- 1_2 and 2_3 x {batch_in_limma}
- 1_2 and 2_3 x {balanced_ruv, all_datasets_ruv}

**Key results (from article):**
- 148 samples integrated (132 first-trimester + 16 second-trimester)
- softImpute + ComBat recovers 61% more genes (21,411 vs 13,250 at 50% coverage)
- Imputation correlation: 0.89 (2_3 cohort), 0.71 (1_2 cohort)
- 99.6% overlap with baseline DEGs; weighted limma confirmed imputation does not introduce false positives

### Fetal Sex Prediction [DONE]

- massiR package used to predict fetal sex for all 148 samples
- 54 male + 78 female (first trimester), 6 male + 10 female (second trimester)
- Sex-stratified DE analysis performed:
  - 1st vs 2nd trimester (sex-agnostic): 368 DEGs
  - 1st trim male vs female: 168 DEGs (129 autosomal)
  - 2nd trim male vs female: 2 DEGs (0 autosomal)
  - 1st vs 2nd trim males: 282 DEGs
  - 1st vs 2nd trim females: 403 DEGs
- Output: `output/yehor_sashko/pack_yehor_1_2_8ds_with_sex_stratified_2026-05-03/`

### Yehor Preprocessing Comparison [DONE]

Compared Sashko (8 datasets) vs Yehor (8-9 datasets) preprocessing:
- Multiple runs in `output/yehor_sashko/` (7ds, 8ds, 9ds variants)
- PCA comparison: `output/yehor_sashko/pca_4_runs_comparison/`
- Leave-one-out validation: `output/yehor_sashko/validation_leave_one_out/`
- GSE37653 intra-dataset batch correction (India/Singapore): `scripts/per_dataset_scripts/GSE37653/`

### GSE37653 Intra-Dataset QC [IN PROGRESS]

- Chip image pseudo-reproductions built from .pair files (`output/qc/GSE37653/`)
- Quality assessment: 8 good, 7 marginal, 10 bad samples (`chip_images_quality_rated.png`)
- Reference NimbleGen datasets downloaded for comparison (`data/raws/nimblegen_examples/`)
- India/Singapore batch effect analyzed via ComBat with 3 sex estimation approaches
- Array design difference documented: 2007 NDF has control probes in grid pattern vs 2010 NDF at edges only

---

## Phase 3: Quality Control [PARTIALLY DONE]

**Chip-level QC (done):**
- [x] Pseudo-image QC for GSE37653 (NimbleGen) - severe spatial artifacts in ~40% of samples
- [x] Reference comparison with GSE55958 and GSE235926
- [x] Per-sample quality rating (good/marginal/bad)

**Study-level QC (not started):**
- [ ] MetaQC package for study quality assessment (IQC, EQC, CQCg/AQCg, CQCp/AQCp)
- [ ] Define exclusion criteria (bottom quartile on multiple metrics)
- [ ] Decide whether to exclude bad GSE37653 samples or keep with batch correction

---

## Phase 4: Validation Framework [PARTIALLY DONE]

**Leave-one-dataset-out cross-validation:**
- [x] Implemented and run (`output/yehor_sashko/validation_leave_one_out/`)
- [x] Also tested without Mikheev dataset (`validation_leave_one_out_mikheev/`)

**Biological validation:**
- [x] Compare with dissertation reference data (`scripts/validation/compare_all_to_disser.R`)
- [ ] Compare with van Uitert 388-gene placenta signature
- [ ] Verify HIF1A pathway enrichment
- [ ] Cross-reference with known placental biology

---

## Phase 5: Baseline Comparison [DONE]

- Dissertation pipeline retained as baseline
- Comparison scripts in `scripts/validation/compare_baseline_to_disser.R`
- 45% DEG overlap between Phase 2B and dissertation (same 6 datasets)
- Differences attributed to sample filtering, ComBat parameters, pipeline implementation

---

## Open Questions

1. **GSE37653 sample quality**: Should bad samples (10 of 25) be excluded? How does this affect batch balance?
2. **Soncin dataset**: Included in 9ds runs but excluded from 8ds — impact on results?
3. **Sex-stratified analysis thresholds**: |logFC| >= log2(1.2) for within-trimester sex comparisons vs |logFC| >= 1 for trimester comparisons — is the lower threshold justified?
4. **Weighted limma**: Should imputed values always be downweighted to 0, or is there a middle ground?

## Key References

**Imputation article & conference abstract:**
- Poliakov, Lykhenko, Obolenska (2026) - Improving integrative analysis by missing value imputation and fetal sex prediction

**Meta-analysis:**
- van Uitert et al. (PLOS ONE, 2015) - Placenta meta-analysis gold standard
- Choi et al. (Bioinformatics, 2003) - Effect-size meta-analysis for genomics
- Breitling et al. (FEBS Letters, 2004) - RankProd method

**Direct merge / Imputation:**
- Mancuso et al. (NAR, 2020) - SampleLASSO cross-platform imputation
- Shabalin et al. (Bioinformatics, 2008) - XPN cross-platform normalization
- Benito et al. (Bioinformatics, 2004) - DWD normalization

## Key Output Locations

| What | Path |
|------|------|
| Runs summary (all methods) | `articles/imputation_article/phase2b_runs_summary.xlsx` |
| Imputation article | `articles/imputation_article/main.pdf` |
| Yehor conference abstract | `articles/yehor_conference_2026/abstract.md` |
| Sex-stratified results | `output/yehor_sashko/pack_yehor_1_2_8ds_with_sex_stratified_2026-05-03/` |
| GSE37653 QC images | `output/qc/GSE37653/chip_images_quality_rated.png` |
| NimbleGen references | `data/raws/nimblegen_examples/` |
| Phase 2B results doc | `scripts/integrative_analysis/research/phase2b_direct_merge_results.md` |
| Leave-one-out validation | `output/yehor_sashko/validation_leave_one_out/` |
