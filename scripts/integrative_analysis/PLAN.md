# Integrative Analysis Improvement Plan

## Problem Statement

Current pipeline loses ~48% of genes through INNER JOIN when merging 8 datasets.
This is unacceptable information loss that likely excludes biologically important genes.

## Recommended Solution: Meta-Analysis Over Direct Merging

### Phase 1: Implement Tiered Gene Coverage Analysis

**Goal:** Replace binary "all or nothing" gene intersection with flexible thresholds.

| Tier | Threshold | Use Case |
|------|-----------|----------|
| Tier 1 | 100% (8/8 studies) | High-confidence biomarker discovery |
| Tier 2 | ≥75% (6/8 studies) | Balanced primary analysis |
| Tier 3 | ≥50% (4/8 studies) | Exploratory analysis |

**Tasks:**
- [ ] Add gene coverage reporting to current pipeline
- [ ] Implement threshold parameter in merging functions
- [ ] Compare DEG results across tiers

### Phase 2: Add Meta-Analysis Methods

**Primary packages to implement:**

1. **DExMA** (Recommended - handles missing genes with imputation)
   - Effect-size meta-analysis
   - Built-in gene imputation via weighted k-nearest samples
   - Flexible gene presence thresholds

2. **RankProd** (Non-parametric alternative)
   - Native missing value handling (excludes from rankings)
   - Platform-independent (uses ranks, not values)
   - Robust to outliers

**Tasks:**
- [ ] Install and test DExMA package
- [ ] Install and test RankProd package
- [ ] Create wrapper functions matching current pipeline interface
- [ ] Compare results with current ComBat+limma approach

### Phase 2B: Direct Merge Improvements (Alternative Path)

**Goal:** If direct merging is preferred over meta-analysis, improve gene recovery through imputation and better normalization.

#### Cross-Platform Gene Imputation

1. **SampleLASSO** (State-of-the-art, Mancuso et al. NAR 2020)
   - Builds sample-specific sparse regression models
   - Outperformed GeneLASSO 91% of the time for cross-platform tasks
   - Biologically interpretable (upweights same tissue samples)
   - GitHub: krishnanlab/Expresto

2. **AffyImpute** (Zhou et al. Bioinformatics 2017)
   - LASSO regression trained on 77,000+ GEO samples
   - Predicts ~10,000 genes unique to HG-U133 Plus 2.0
   - Spearman correlation 0.90-0.96 between imputed and measured

3. **Matrix completion** (softImpute - Hastie et al. JMLR 2015)
   - Nuclear norm minimization
   - Exploits low-rank structure of expression matrices

#### Alternative Normalization Methods (Instead of ComBat)

| Method | Description | Best For |
|--------|-------------|----------|
| **XPN** | Cross-platform normalization via gene/sample clustering | Equal-sized treatment groups |
| **DWD** | Distance Weighted Discrimination, margin-based | Unequal group sizes, robust |
| **CuBlock** | Cubic polynomial per probe cluster | Sample-by-sample processing |
| **crossmeta** | Effect-size combination (not true merging) | Different gene universes |

**Tasks:**
- [ ] Evaluate SampleLASSO for Affymetrix-Illumina imputation
- [x] Test DWD as ComBat alternative (result: too aggressive, not recommended)
- [x] Benchmark imputation quality at 12.3% missingness (correlation: 0.89)
- [x] Document validation: imputed vs measured correlation (see research/phase2b_direct_merge_results.md)
- [x] Implement softImpute matrix completion

**Results (2026-01-18):**
- softImpute + ComBat achieves 0.89 correlation and recovers 61% more genes
- DWD normalization too aggressive (eliminates biological signal)
- See `research/phase2b_direct_merge_results.md` for full analysis
- **Permanent data folder**: `output/phase2b_term_vs_2trim_final/`

**Caveats:**
- At 48% missingness, beyond typical validation scenarios (performance degrades >20%)
- Imputation validity depends heavily on validation
- Meta-analysis often safer than forcing integration through imputation

---

### Phase 3: Quality Control Integration

**Implement MetaQC package for study quality assessment:**

- IQC: Internal coexpression homogeneity
- EQC: External pathway consistency
- CQCg/AQCg: DE gene accuracy/consistency
- CQCp/AQCp: Pathway accuracy/consistency

**Tasks:**
- [ ] Add MetaQC to preprocessing step
- [ ] Define exclusion criteria (bottom quartile on multiple metrics)
- [ ] Document all exclusion reasons

### Phase 4: Validation Framework

**Leave-one-dataset-out cross-validation:**
- Train on 7 datasets, validate on held-out study
- Repeat for all 8 combinations
- More realistic than within-study CV

**Biological validation:**
- [ ] Compare with van Uitert 388-gene placenta signature
- [ ] Verify HIF1A pathway enrichment
- [ ] Cross-reference with known placental biology

### Phase 5: Keep Current Pipeline as Baseline

Retain `integrative_analysis_disser_pipeline/` for:
- Comparison with new methods
- Backward compatibility
- Published dissertation methodology

---

## Package Requirements

```r
# Meta-analysis packages (Phase 2)
BiocManager::install(c(
  "DExMA",      # Meta-analysis with imputation
  "RankProd",   # Rank-based meta-analysis
  "crossmeta",  # Automated multi-platform meta-analysis
  "MetaDE",     # Multiple meta-analysis methods
  "metafor"     # General meta-analysis, forest plots
))

# Direct merge improvements (Phase 2B)
# SampleLASSO: devtools::install_github("krishnanlab/Expresto")
# CONOR (XPN, DWD): install.packages("CONOR")
install.packages(c(
  "softImpute", # Matrix completion
  "glmnet"      # LASSO for imputation
))

# Quality control (Phase 3)
BiocManager::install("MetaQC")
```

## Directory Structure (Proposed)

```
scripts/integrative_analysis/
├── PLAN.md                              # This file
├── research/                            # Background research
├── integrative_analysis_disser_pipeline/  # Current approach (baseline)
├── meta_analysis_pipeline/              # Phase 2: Meta-analysis methods
│   ├── 01_quality_control.R             # MetaQC integration
│   ├── 02_within_study_de.R             # Per-study limma
│   ├── 03_effect_size_meta.R            # DExMA wrapper
│   ├── 04_rank_meta.R                   # RankProd wrapper
│   └── 05_validation.R                  # Cross-study validation
├── direct_merge_pipeline/               # Phase 2B: Improved direct merge
│   ├── 01_imputation.R                  # SampleLASSO, softImpute
│   ├── 02_cross_platform_norm.R         # XPN, DWD alternatives
│   └── 03_merge_with_imputed.R          # Merge after imputation
└── comparison/                          # Compare all approaches
    └── compare_methods.R
```

## Key References

**Meta-analysis:**
- van Uitert et al. (PLOS ONE, 2015) - Placenta meta-analysis gold standard
- Choi et al. (Bioinformatics, 2003) - Effect-size meta-analysis for genomics
- Breitling et al. (FEBS Letters, 2004) - RankProd method
- Villatoro-García et al. (Mathematics, 2022) - DExMA package

**Direct merge / Imputation:**
- Mancuso et al. (NAR, 2020) - SampleLASSO cross-platform imputation
- Zhou et al. (Bioinformatics, 2017) - AffyImpute
- Shabalin et al. (Bioinformatics, 2008) - XPN cross-platform normalization
- Benito et al. (Bioinformatics, 2004) - DWD normalization
- Rudy & Valafar (BMC Bioinformatics, 2011) - CONOR benchmark

## Success Metrics

1. **Gene coverage:** Retain >80% of genes (vs current 52%)
2. **DEG overlap:** >70% overlap with current method at Tier 1
3. **Validation:** Results replicate in leave-one-out analysis
4. **Biology:** Enrichment for known placental pathways

## Test-Driven Development (TDD) Approach

All new pipeline components should follow TDD principles:

### Testing Strategy

1. **Write tests BEFORE implementation**
   - Define expected behavior first
   - Tests document requirements
   - Prevents scope creep

2. **Test categories:**

| Category | What to Test | Example |
|----------|--------------|---------|
| Unit tests | Individual functions | `test_gene_coverage_calculation()` |
| Integration tests | Module interactions | `test_tiered_merge_with_combat()` |
| Regression tests | Comparison with baseline | `test_deg_overlap_with_disser()` |
| Biological validation | Known biology | `test_hif1a_pathway_enrichment()` |

### Test Files Structure

```
tests/
├── test_phase1_tiered_coverage.R
│   ├── test_gene_presence_matrix()
│   ├── test_coverage_at_thresholds()
│   └── test_tier_filtering()
├── test_phase2_meta_analysis.R
│   ├── test_dexma_wrapper()
│   ├── test_rankprod_wrapper()
│   └── test_meta_vs_baseline_overlap()
├── test_phase2b_imputation.R
│   ├── test_imputation_correlation()
│   └── test_imputed_vs_measured()
├── test_phase3_qc.R
│   └── test_metaqc_exclusion_criteria()
└── test_validation.R
    ├── test_leave_one_out_cv()
    └── test_van_uitert_signature_overlap()
```

### Example Test (Phase 1)

```r
# tests/test_phase1_tiered_coverage.R

test_that("gene_coverage returns correct percentages", {
  # Setup: 3 datasets with overlapping genes
  dataset1 <- c("A", "B", "C", "D")
  dataset2 <- c("A", "B", "E")
  dataset3 <- c("A", "C", "E", "F")

  coverage <- calculate_gene_coverage(list(dataset1, dataset2, dataset3))

  expect_equal(coverage["A"], 1.0)      # In all 3

expect_equal(coverage["B"], 2/3)      # In 2 of 3
  expect_equal(coverage["F"], 1/3)      # In 1 of 3
})

test_that("tier filtering returns correct gene counts", {
  coverage <- c(A=1.0, B=0.75, C=0.5, D=0.25)

  expect_equal(length(filter_by_tier(coverage, 1.0)), 1)   # Only A
  expect_equal(length(filter_by_tier(coverage, 0.75)), 2)  # A, B
  expect_equal(length(filter_by_tier(coverage, 0.5)), 3)   # A, B, C
})
```

### TDD Workflow

1. **Red**: Write failing test for new feature
2. **Green**: Implement minimum code to pass
3. **Refactor**: Clean up while keeping tests green
4. **Repeat**: Next feature

### Continuous Integration

- Run `testthat::test_dir("tests/")` before committing
- All tests must pass before merging new pipeline code
- Track test coverage with `covr` package

---

## Timeline

1. Phase 1 (Tiered coverage): Foundation for flexibility
2. Phase 2 (Meta-analysis): Core improvement
3. Phase 3 (QC): Study quality filtering
4. Phase 4 (Validation): Ensure robustness
5. Phase 5 (Comparison): Document improvements
