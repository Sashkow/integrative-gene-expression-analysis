## Context

The pipeline currently selects RUV negative control genes via a top-down approach (`build_ruv_control_genes.R`): start from the Eisenberg-Levanon housekeeping list, then filter by empirical DE status. This change adds a complementary bottom-up approach: select genes empirically from Mikheev (GSE9984) and Soncin (GSE100051) — the two datasets that span both 1st and 2nd trimesters — based on expression stability, then cross-validate against literature.

Both datasets are already preprocessed as ENTREZID-keyed, log2-transformed, protein-coding TSVs in `data/mapped/yehor/8_yehor_preprocessed_datasets_enriched_sashko/`.

The existing `build_ruv_control_genes.R` and its literature verdicts (lines 68-99) provide a reference for cross-checking results.

## Goals / Non-Goals

**Goals:**
- Select empirical RUV negative control gene candidates from expression data
- Apply the 25th-90th percentile expression filter to avoid noise floor and probe saturation
- Rank by coefficient of variation (CV) within the expression band
- Exclude genes that are differentially expressed in 1_2 contrast
- Cross-reference with Eisenberg-Levanon housekeeping list
- Cross-reference with literature-known 1_2 DEGs (qPCR/blot validated)
- Produce annotated output (per-gene table + plain ENTREZ list)

**Non-Goals:**
- Not replacing the existing `build_ruv_control_genes.R` — this is a complementary approach
- Not running RUV itself — only selecting candidate genes
- Not optimizing RUV k parameter — that's a separate concern
- Not handling the 2_3 contrast — focused on 1_2 only (can be extended later)

## Decisions

### 1. Use both Mikheev and Soncin, not just one

**Decision**: Compute stability metrics independently in each dataset, then take the intersection of candidates.

**Rationale**: A gene that is stable in both an Affymetrix dataset (GSE9984) and an Illumina dataset (GSE100051) is a stronger candidate than one stable in only one platform. This also provides cross-platform validation of stability.

**Alternative considered**: Pool all samples into one matrix before computing CV. Rejected because batch effects between datasets would inflate variance of stable genes and deflate variance of batch-correlated genes.

### 2. Use CV (coefficient of variation) not raw variance

**Decision**: Rank by CV = SD / mean within the expression band, not by raw SD or variance.

**Rationale**: Raw variance has a mean-variance relationship in log-expression data — genes with higher average expression tend to have slightly higher variance. CV normalizes for this. Since we already filter by expression percentile, this is a secondary adjustment but makes ranking fairer within the band.

### 3. Expression filter: 25th-90th percentile

**Decision**: Filter to genes with average expression between the 25th and 90th percentile, computed per-dataset.

**Rationale**: Below the 25th percentile, signal approaches the noise floor and apparent low variance may reflect measurement insensitivity rather than biological stability. Above the 90th percentile, probe saturation compresses fold changes and may create artificial stability. The 25th-90th band retains ~65% of genes.

### 4. DE exclusion source

**Decision**: Use the pipeline's existing DE results from the full 6-dataset softImpute + ComBat-ref run as the primary DE exclusion source, supplemented by per-dataset DE tests within Mikheev and Soncin.

**Rationale**: The full pipeline DE results (`difexp_softimpute_combat_ref.tsv`) have the most statistical power (117 samples). Per-dataset DE adds dataset-specific sensitivity. A gene excluded by any source is excluded.

### 5. Literature DEG list

**Decision**: Curate a list of genes known from qPCR/Western blot studies to be differentially expressed in 1st-vs-2nd trimester human placenta. Source from Prater 2021 top DEGs (validated by qPCR in their study) and other published placental DE studies.

**Rationale**: Empirical filters might miss genes whose DE is real but too small to reach significance in our pipeline. Literature-validated DEGs should be excluded regardless of our pipeline's verdict.

## Risks / Trade-offs

- **Too few candidates**: If the intersection of all filters is very small, the gene list may be insufficient for RUV. Mitigation: report counts at each filter step so the user can relax thresholds.
- **Platform-specific stability**: A gene stable on Affymetrix might not be stable on Illumina due to probe design differences. Mitigation: requiring stability in both datasets handles this, at the cost of a smaller final list.
- **Literature list incompleteness**: We may not find an exhaustive list of literature-validated 1_2 DEGs. Mitigation: this filter is additive (it only excludes genes, never includes); missing a few literature DEGs is acceptable since the empirical filters catch most of them.
