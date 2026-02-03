# Meta-Analysis vs Data Merging: Comparative Analysis

**Date**: 2026-01-17
**Comparison**: Term vs Second Trimester placental gene expression

## Overview

This document compares two approaches for integrating gene expression data from multiple GEO studies:

1. **Data Merging (Baseline)**: ComBat batch correction followed by limma differential expression
2. **Meta-Analysis**: Per-study analysis followed by effect size combination (DExMA, RankProd, Metafor)

---

## Datasets

### Studies Available

| Study | Platform | Second Trimester | Term | Has Both Conditions |
|-------|----------|-----------------|------|---------------------|
| GSE100051 | Affymetrix | 7 | 5 | Yes |
| GSE9984 | Affymetrix | 4 | 4 | Yes |
| GSE22490 | Affymetrix | 1 | 0 | No |
| GSE37901 | Affymetrix | 4 | 0 | No |
| GSE35574 | Affymetrix | 0 | 26 | No |
| GSE55439 | Illumina | 0 | 15 | No |
| GSE6573 | Affymetrix | 0 | 1 | No |
| GSE73374 | Affymetrix | 0 | 12 | No |
| GSE73685 | Affymetrix | 0 | 12 | No |
| **Total** | | **16** | **75** | |

### Data Used by Each Method

| Method | Studies Used | Samples | Rationale |
|--------|-------------|---------|-----------|
| **Baseline (ComBat+limma)** | All 9 | 91 (16 vs 75) | Merges all data after batch correction |
| **Metafor** | 2 studies* | 20 (11 vs 9) | Requires samples in both conditions for effect sizes |
| **DExMA** | 2 studies | 20 (11 vs 9) | Requires samples in both conditions |
| **RankProd** | 2 studies | 20 (11 vs 9) | Requires samples in both conditions |

**Key limitation**: Only GSE100051 and GSE9984 have samples in both gestational age groups. All three meta-analysis methods effectively use only these 2 studies.

*Metafor attempts to run on all 7 loaded studies, but unbalanced studies produce "Coefficients not estimable" warnings and contribute no effect sizes.

---

## Methods Description

### Baseline: ComBat + limma
- Merges expression matrices from all studies
- Applies ComBat batch correction to remove study effects
- Runs limma differential expression on merged data
- Reports log2 fold change (logFC) and adjusted p-values

### Meta-Analysis Methods

#### Metafor (Effect-Size Random Effects)
- Runs limma on each study separately to compute effect sizes
- Combines effect sizes using random effects model (DerSimonian-Laird)
- **Requires balanced studies**: Unbalanced studies produce "Coefficients not estimable" and contribute nothing
- In practice, only GSE100051 and GSE9984 contributed effect sizes
- Reports combined logFC and FDR

#### DExMA (Differential Expression Meta-Analysis)
- Requires balanced studies (samples in both conditions)
- Uses Standardized Mean Difference (SMD) as effect size
- Applies random effects model with optional imputation
- Reports SMD (labeled as "logFC" in output) and FDR

#### RankProd (Rank Product)
- Non-parametric rank-based method
- Requires balanced studies
- Identifies consistently up/down-regulated genes across studies
- Reports Percentage of False Positives (PFP)

---

## Important: SMD vs logFC

**Critical finding**: DExMA reports Standardized Mean Difference (SMD), not log fold change.

```
logFC = mean(case) - mean(control)           # Units: log2 expression
SMD   = (mean(case) - mean(control)) / SD    # Units: standard deviations (Cohen's d)
```

| logFC | SD | SMD | Interpretation |
|-------|-----|-----|----------------|
| 0.5 | 0.3 | 1.67 | Small fold change, large effect (low variance) |
| 2.0 | 4.0 | 0.5 | Large fold change, small effect (high variance) |
| 1.0 | 1.0 | 1.0 | Equal |

**Implication**: A gene with modest logFC but consistent expression (low variance) will have large SMD. This explains why all DExMA significant genes have |SMD| > 1.

---

## Results Summary

### Raw Counts (Before Fair Comparison)

| Method | Total Genes | Significant (FDR < 0.05) | With |logFC| > 1 |
|--------|-------------|-------------------------|-----------------|
| Baseline | 13,250 | 7,989 | 534 |
| Metafor | 18,039 | 4,148 | 417 |
| DExMA | 16,412 | 3,177 | 3,177* |
| RankProd | 16,412 | 2,292 | N/A |
| High Confidence | 1,491 | 1,491 | 1,491* |

*All DExMA significant genes have |SMD| > 1 because SMD tends to be larger than logFC

### Fair Comparison (Using SMD for Both)

To enable fair comparison, baseline logFC was converted to SMD using the t-statistic:
```
SMD = t × sqrt(1/n1 + 1/n2)
```

Where:
- **t** = t-statistic from limma output (column "t" in difexp_all.csv) - measures significance of difference
- **n1** = number of Term samples = 75
- **n2** = number of Second Trimester samples = 16

This formula derives from: t = (mean1 - mean2) / SE, where SE (Standard Error) = pooled_SD × sqrt(1/n1 + 1/n2), therefore SMD = logFC / pooled_SD = t × sqrt(1/n1 + 1/n2).

**With FDR < 0.05 and |SMD| > 1:**

| Metric | Baseline | DExMA |
|--------|----------|-------|
| Significant genes | 5,287 | 3,016 |
| Overlap | 2,775 | 2,775 |
| Recovery rate | 52.5% | 92.0% |

Recovery rate calculation:
- Baseline: 2,775 / 5,287 = 52.5% (% of baseline genes also found in DExMA)
- DExMA: 2,775 / 3,016 = 92.0% (% of DExMA genes also found in baseline)

**SMD correlation between methods: r = 0.941**

### Comparison Using logFC (Converting DExMA SMD to logFC)

To compare apples-to-apples, DExMA's SMD was converted back to logFC using:
```
logFC_dexma = SMD_dexma × pooled_SD
```
Where pooled_SD is calculated from baseline: pooled_SD = logFC_baseline / SMD_baseline

**logFC correlation (baseline vs converted DExMA): r = 0.929**

#### All Genes (FDR < 0.05 only)

| Metric | Baseline | DExMA | Overlap |
|--------|----------|-------|---------|
| Significant genes | 7,901 | 3,016 | 2,988 |
| Recovery rate | - | - | 37.8% of baseline |

#### With |logFC| > 1 Filter

| Metric | Baseline | DExMA (converted) | Overlap |
|--------|----------|-------------------|---------|
| FDR<0.05 + \|logFC\|>1 | 532 | 369 | 243 |
| Recovery rate | - | - | 45.7% of baseline |

The 243 genes significant in both methods with |logFC| > 1 in both represent the most robust large-effect findings.

---

## Key Findings

### 1. High Agreement Between Methods
- Effect size correlation: **r = 0.941** (SMD) and **r = 0.987** (logFC with Metafor)
- Direction agreement: **100%** for all overlapping significant genes
- The methods are finding the same biological signal

### 2. DExMA is a Conservative Subset
- 92% of DExMA significant genes are also significant in baseline
- DExMA finds fewer genes because it only uses 2 balanced studies
- The genes it finds are highly reliable (cross-study validated)

### 3. Baseline Has More Power but Less Validation
- Uses all 91 samples → more statistical power
- Finds more significant genes (5,287 vs 3,016 at |SMD| > 1)
- But findings are not cross-validated across independent studies
- Unbalanced comparison (75 Term vs 16 Second Trimester)

### 4. High Confidence Gene Set
- 1,491 genes significant in both DExMA AND RankProd
- 100% direction agreement between methods
- 1,450 of these (97%) are also significant in baseline
- **Most robust finding set for downstream analysis**

### 5. Unique Discoveries
- 41 genes significant in high-confidence meta-analysis but NOT in baseline
- These may represent effects that are consistent across studies but masked in merged analysis
- 4,322 genes significant in baseline but NOT in Metafor
- These may be study-specific effects or genes with inconsistent patterns

**41 genes unique to meta-analysis:**
```
ERVMER34-1
MYZAP
C1D
USP43
SPAAR
DEFA4
TMEM273
CALHM4
ZBTB38
TMEM256
STEAP1B
IAH1
GPSM2
RHOD
H2AX
HBM
HMGB1
HSPA6
CENPS
SKIDA1
SMIM1
C20orf202
NICOL1
TOMM5
CKLF
TMEM132A
TMIGD3
BGN
CCL4
SLC18A3
TAF9
XG
PXDN
TAF1D
PGBD5
SVEP1
H4C8
PIR
PGLYRP1
NREP
EEF1E1
```
(no significan enrichment detected with stringdb)
---

## Recommendations

### For Discovery (Exploratory Analysis)
Use **baseline (ComBat+limma)** results:
- Maximum statistical power
- Largest gene set for pathway/enrichment analysis
- 5,287 genes at FDR < 0.05, |SMD| > 1

### For Validation (High Confidence)
Use **meta-analysis high confidence** set:
- 1,491 genes significant in both DExMA and RankProd
- Cross-validated across independent studies
- Best candidates for experimental follow-up

### For Publication
Report both:
1. Baseline analysis as primary (more power)
2. Meta-analysis as validation (independent replication)
3. Highlight the 1,450 genes significant in both approaches

---

## Output Files

```
output/phase2_meta/
├── metafor_results.tsv          # Metafor full results
├── metafor_significant.tsv      # FDR < 0.05
├── dexma_results.tsv            # DExMA full results (SMD as effect size)
├── dexma_significant.tsv        # FDR < 0.05
├── rankprod_results.tsv         # RankProd full results
├── rankprod_significant.tsv     # PFP < 0.05
├── high_confidence_degs.tsv     # Significant in both DExMA + RankProd
├── combined_results.tsv         # All methods merged
├── comparison/
│   ├── baseline_with_smd.csv    # Baseline with calculated SMD
│   ├── merged_smd_comparison.csv # Side-by-side SMD comparison
│   ├── comparison_summary.csv   # Summary statistics
│   ├── comparison_report.txt    # Text report
│   └── *.pdf                    # Correlation plots, Venn diagrams
└── summary.txt                  # Phase 2 execution summary
```

---

## Technical Notes

### Why DExMA/RankProd Use Fewer Studies
These methods require each study to have samples in both conditions to calculate within-study effect sizes. In this dataset, only 2 of 9 studies (GSE100051, GSE9984) meet this criterion.

### Metafor Gene Coverage
Metafor reports more genes (18,039) than DExMA/RankProd (16,412) because it uses a "min 4 studies" threshold for gene inclusion. However, since only 2 studies contribute valid effect sizes (the balanced ones), the additional genes come from having measurements in multiple unbalanced studies - which don't actually contribute to the meta-analysis. The effective analysis is based on the same 2 balanced studies as DExMA/RankProd.

### SMD Calculation from Baseline
The baseline reports t-statistics from limma. SMD was calculated as:
```r
SMD = t × sqrt(1/n_term + 1/n_second)
# where n_term = 75, n_second = 16
```

This approximation is valid for two-group comparisons and allows direct comparison with DExMA's SMD values.
