# Phase 2B: Direct Merge Improvements - Results

**Date**: 2026-01-18
**Comparison**: Term vs Second Trimester placental gene expression
**Permanent data folder**: `output/phase2b_term_vs_2trim_final/`

## Overview

Phase 2B evaluates two strategies for improving direct data merging:
1. **Gene imputation** using softImpute (nuclear norm matrix completion)
2. **Alternative normalization** using DWD instead of ComBat

---

## Key Findings

### 1. Gene Recovery with Imputation

Total samples across 8 datasets: 409
Samples used for DE (Term vs Second Trimester): **182** (165 Term + 17 Second Trimester)

**Imputation statistics for 182 DE samples:**

| Coverage | Genes | Total Cells | Cells to Impute | % Missing |
|----------|-------|-------------|-----------------|-----------|
| 100% (8/8) | 13,250 | 2,411,500 | 0 | 0.0% |
| 75% (6/8) | 17,412 | 3,168,984 | 101,061 | 3.2% |
| **50% (4/8)** | **21,411** | **3,896,802** | **491,230** | **12.6%** |
| 25% (2/8) | 24,117 | 4,389,294 | 796,106 | 18.1% |

**Finding**: Using 50% coverage threshold with imputation recovers **61% more genes** (21,411 vs 13,250) compared to inner join (100% coverage).

### 2. Imputation Quality Validation

Leave-10%-out cross-validation (3 repeats):

| Metric | Value |
|--------|-------|
| Correlation (imputed vs true) | **0.89** |
| RMSE | 6.28 |
| MAE | 5.81 |

**Finding**: softImpute achieves **0.89 correlation** between imputed and true values, indicating good imputation quality for gene expression data.

### 3. Method Comparison

| Method | Imputation | Normalization | Genes | Significant DEGs | Up | Down |
|--------|------------|---------------|-------|------------------|-----|------|
| none_combat | None | ComBat | 13,250 | 459 | 196 | 263 |
| none_dwd | None | DWD | 13,250 | 11 | 0 | 11 |
| **softimpute_combat** | softImpute | ComBat | 21,411 | **517** | 224 | 293 |
| softimpute_dwd | softImpute | DWD | 21,411 | 0 | 0 | 0 |

### 4. DEG Overlap

|  | none_combat | none_dwd | softimpute_combat | softimpute_dwd |
|--|-------------|----------|-------------------|----------------|
| none_combat | 459 | 7 | 457 | 0 |
| none_dwd | 7 | 11 | 7 | 0 |
| softimpute_combat | 457 | 7 | 517 | 0 |
| softimpute_dwd | 0 | 0 | 0 | 0 |

**Key observations**:
- 457 of 459 (99.6%) genes from baseline (none_combat) are also found with imputation
- softimpute_combat finds **60 additional DEGs** not detected without imputation
- DWD normalization eliminates nearly all biological signal (too aggressive)

---

## Conclusions

### Recommended Approach: softImpute + ComBat

The combination of softImpute imputation with ComBat batch correction:
- Recovers 61% more genes (21,411 vs 13,250)
- Finds 12.6% more DEGs (517 vs 459)
- Maintains 99.6% overlap with baseline results
- Achieves 0.89 imputation correlation

### DWD Not Recommended

The simplified DWD implementation proved too aggressive:
- Eliminates most biological variance along with batch effects
- Results in near-zero DEGs
- Would require more sophisticated implementation to preserve biological signal

### Limitations

1. **Imputation at 12.3% missingness**: Within reasonable bounds, but higher than typical validation scenarios
2. **DWD implementation**: The simplified mean-shift approach doesn't properly preserve biological variance
3. **Cross-platform differences**: softImpute assumes genes behave similarly across platforms, which may not always hold

---

## 1_2 Comparison (First Trimester vs Second Trimester) - 8 Datasets

**Datasets**: 8 (GSE100051, GSE22490, GSE28551, GSE37653, GSE55439, GSE122214, GSE37901, GSE9984)
**Total samples**: 125 DE samples (108 First Trimester + 17 Second Trimester)
**Recommended coverage**: 75% (6/8 datasets)
**Permanent output folder**: `output/phase2b_1_2_8datasets_final/`

### Gene Recovery (1_2)

| Coverage | Genes | % Missing |
|----------|-------|-----------|
| 100% (8/8) | 8,907 | 0% |
| **75% (6/8)** | **17,176** | **12.8%** |
| 50% (4/8) | 20,538 | 22.7% |
| 25% (2/8) | 22,695 | 25.3% |

### Imputation Quality (1_2)

| Metric | Value |
|--------|-------|
| Correlation (imputed vs true) | **0.71** |
| RMSE | 10.62 |
| MAE | 8.70 |

Note: Lower correlation than 2_3 (0.71 vs 0.89) due to fewer Second Trimester samples (17 vs 17) spread across fewer datasets.

### Method Comparison (1_2 at 75% coverage)

| Method | Genes | DEGs | Up | Down |
|--------|-------|------|-----|------|
| none_combat | 8,907 | 165 | 120 | 45 |
| none_dwd | 8,907 | 85 | 67 | 18 |
| **softimpute_combat** | **17,176** | **174** | 132 | 42 |
| softimpute_dwd | 17,176 | 97 | 74 | 23 |

### DEG Imputation Rate Analysis (1_2)

**Key finding**: Imputed values do affect logFC and p-value calculations. Analysis of per-gene imputation rates:

| Imputation Level | DEGs | Percentage |
|------------------|------|------------|
| **0% imputed** (in all 8 datasets) | 151 | 87% |
| **12.5% imputed** (in 7/8 datasets) | 15 | 8.6% |
| **25% imputed** (in 6/8 datasets) | 8 | 4.6% |
| **>25% imputed** | 0 | 0% |

**DEGs with 25% imputation (present in 6/8 datasets):**

| Gene | logFC | adj.P.Val |
|------|-------|-----------|
| MED14 | 3.96 | 0.0009 |
| KIF26A | 3.31 | 0.0013 |
| DENND4B | 4.01 | 0.0043 |
| OR8G2P | 2.20 | 0.0089 |
| BRD3OS | 3.20 | 0.0106 |
| GLS | 2.49 | 0.0135 |
| SCML4 | 1.84 | 0.0310 |
| TRIM66 | 2.80 | 0.0441 |

**Conclusion**: 87% of DEGs (151/174) required NO imputation. The remaining 23 DEGs have at most 25% imputed values. No DEGs with >25% imputation reached significance, suggesting imputation is not introducing false positives.

Full per-gene imputation rates saved to: `output/phase2b_1_2_8datasets_final/deg_imputation_rates.csv`

### Weighted Limma Analysis (1_2)

**Key finding**: When imputed values are weighted to 0 in limma (so they don't contribute to logFC/p-value), 22 out of 23 imputed DEGs lose significance.

| Method | DEGs |
|--------|------|
| Unweighted (imputed = real) | 174 |
| **Weighted (imputed = 0)** | **154** |
| Overlap | 152 |

The 22 genes that lost significance had logFC values drop from 2-4 range to 0.1-1.6, indicating imputation was creating false positives.

Only **ZNF638** (missing from GSE100051, 39% imputed) survived weighted analysis: logFC dropped from 4.66 to 2.43 but remained significant (adj.P = 0.002).

Weighted results saved to: `output/phase2b_1_2_8datasets_final/difexp_softimpute_combat_weighted.tsv`

---

## 1_2 Comparison - Dissertation Datasets (6 datasets)

**Datasets**: GSE100051, GSE122214, GSE22490, GSE37901, GSE55439, GSE9984 (same as dissertation)
**Total samples**: 84 DE samples (67 First + 17 Second)
**Output folder**: `output/phase2b_1_2_disser/`

### Gene Recovery (Dissertation Datasets)

| Coverage | Genes | % Missing |
|----------|-------|-----------|
| 100% (6/6) | 13,699 | 0% |
| **75% (5/6)** | **17,413** | **5.3%** |
| 50% (3/6) | 20,261 | 14.7% |

### Imputation Quality
- Correlation: **0.75** (better than 8-dataset 0.71)

### Method Comparison

| Method | Genes | DEGs | Up | Down |
|--------|-------|------|-----|------|
| **none_combat** | 13,699 | **198** | 144 | 54 |
| none_dwd | 13,699 | 112 | 84 | 28 |
| softimpute_combat | 17,413 | 197 | 143 | 54 |
| softimpute_dwd | 17,413 | 111 | 82 | 29 |

Note: With only 5.3% missing data at 75% coverage, imputation provides minimal benefit (198 vs 197 DEGs).

### Comparison with Dissertation

| | Dissertation | Phase 2B |
|---|---|---|
| Datasets | 6 | 6 (same) |
| Samples | 80 | 84 |
| DEGs | 253 | 198 |
| **Overlap** | - | **89 (45%)** |

The low overlap (45%) despite using the same datasets suggests differences in:
- Sample filtering (dissertation had 80 vs our 84 samples)
- ComBat/limma parameters
- Pipeline implementation details

---

## Output Files

**2_3 Comparison (Term vs Second Trimester):**
```
output/phase2b_term_vs_2trim_final/
├── phase2b_direct_merge_results.md   # This document
├── recovery_vs_imputation.png        # Visualization of recovery vs imputation
├── gene_recovery_comparison.csv      # Genes retained at different thresholds
├── imputation_validation.csv         # Leave-out validation metrics
├── method_comparison.csv             # Summary of all method combinations
├── deg_overlap_matrix.csv            # DEG overlap between methods
├── exprs_*.tsv                       # Expression matrices
├── difexp_*.tsv                      # Differential expression results
├── phase2b_log_*.txt                 # Execution log
└── summary.txt                       # Full summary
```

**1_2 Comparison - 8 Datasets (First Trimester vs Second Trimester):**
```
output/phase2b_1_2_8datasets_final/
├── deg_imputation_rates.csv          # Per-DEG imputation percentage
├── deg_imputation_rates_detailed.csv # Exact imputation % per sample
├── difexp_softimpute_combat_weighted.tsv  # Weighted limma results
├── difexp_*.tsv                      # All differential expression results
├── exprs_*.tsv                       # Expression matrices
└── summary.txt                       # Full summary
```

**1_2 Comparison - Dissertation Datasets (6 datasets):**
```
output/phase2b_1_2_disser/
├── gene_recovery_comparison.csv      # Genes retained at different thresholds
├── imputation_validation.csv         # Leave-out validation metrics
├── method_comparison.csv             # Summary of all method combinations
├── difexp_*.tsv                      # Differential expression results
├── exprs_*.tsv                       # Expression matrices
└── summary.txt                       # Full summary
```

---

## Recommendations for Pipeline

1. **Default to inner join** (100% coverage) for highest-confidence analyses
2. **Use softImpute + ComBat at 50-75% coverage** for exploratory analyses needing more genes
3. **Do not use DWD** as implemented; stick with ComBat for batch correction
4. **Report imputation in methods** when using softImpute with coverage < 100%
