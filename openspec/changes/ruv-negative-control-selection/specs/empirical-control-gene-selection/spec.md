## ADDED Requirements

### Requirement: Expression level filtering
The script SHALL load per-dataset expression matrices for GSE9984 (Mikheev) and GSE100051 (Soncin) and compute per-gene average expression across all samples within each dataset. Genes with average expression below the 25th percentile or above the 90th percentile (computed per-dataset) SHALL be excluded from further consideration.

#### Scenario: Gene below noise floor excluded
- **WHEN** a gene has average expression below the 25th percentile in a dataset
- **THEN** the gene is excluded from that dataset's candidate list

#### Scenario: Gene above saturation ceiling excluded
- **WHEN** a gene has average expression above the 90th percentile in a dataset
- **THEN** the gene is excluded from that dataset's candidate list

### Requirement: Variance-based ranking
Within the expression-filtered gene set, the script SHALL compute the coefficient of variation (CV = SD / mean) for each gene across all samples in each dataset. Genes SHALL be ranked by ascending CV (lowest CV = most stable).

#### Scenario: CV computed per dataset
- **WHEN** expression filtering is complete for a dataset
- **THEN** CV is computed for each remaining gene using all samples in that dataset
- **THEN** genes are ranked by ascending CV

### Requirement: Cross-dataset intersection
The script SHALL compute stability metrics independently for each dataset and retain only genes that pass the expression filter in both datasets. The final stability rank SHALL be the average CV rank across both datasets.

#### Scenario: Gene stable in both datasets
- **WHEN** a gene passes the expression filter in both GSE9984 and GSE100051 and has low CV in both
- **THEN** the gene is retained as a candidate

#### Scenario: Gene stable in only one dataset
- **WHEN** a gene passes the expression filter in only one dataset or has high CV in one
- **THEN** the gene is excluded from the final candidate list

### Requirement: Differential expression exclusion
The script SHALL exclude any gene that is differentially expressed (adj.P.Val < 0.05) in the full pipeline DE results (softImpute + ComBat-ref, 1_2 comparison). The script SHALL also perform per-dataset limma DE tests within GSE9984 and GSE100051 (1st vs 2nd trimester) and exclude genes significant at FDR < 0.05 in either dataset.

#### Scenario: Gene is DE in full pipeline
- **WHEN** a gene has adj.P.Val < 0.05 in the full pipeline difexp table
- **THEN** the gene is excluded regardless of its CV rank

#### Scenario: Gene is DE within one dataset
- **WHEN** a gene has FDR < 0.05 in a per-dataset 1T-vs-2T limma test
- **THEN** the gene is excluded regardless of its CV rank

### Requirement: Housekeeping list cross-reference
The script SHALL annotate each candidate gene with whether it appears in the Eisenberg-Levanon (2013) housekeeping gene list. This annotation is informational — it does not filter genes in or out.

#### Scenario: Candidate is in housekeeping list
- **WHEN** a candidate gene's ENTREZ ID maps to a symbol in the Eisenberg-Levanon list
- **THEN** the gene is annotated with housekeeping_flag = TRUE

#### Scenario: Candidate is not in housekeeping list
- **WHEN** a candidate gene is not in the Eisenberg-Levanon list
- **THEN** the gene is annotated with housekeeping_flag = FALSE

### Requirement: Literature DEG exclusion
The script SHALL maintain a curated list of genes known from qPCR, Western blot, or other targeted experiments to be differentially expressed in 1st-vs-2nd trimester human placenta. Any gene in this list SHALL be excluded from the final candidate set.

#### Scenario: Gene is a known literature DEG
- **WHEN** a candidate gene appears in the literature DEG exclusion list
- **THEN** the gene is excluded from the final candidate list

### Requirement: Output files
The script SHALL produce:
1. An annotated per-gene table (TSV or XLSX) with columns: ENTREZID, SYMBOL, mean_expr_GSE9984, mean_expr_GSE100051, CV_GSE9984, CV_GSE100051, avg_CV_rank, pipeline_FDR, per_dataset_DE, housekeeping_flag, literature_DEG, final_decision, exclusion_reason
2. A plain text file of ENTREZ IDs for the final selected genes (one per line), suitable for direct input to the RUV pipeline
3. A summary printed to console with gene counts at each filter step

#### Scenario: All outputs produced
- **WHEN** the script completes successfully
- **THEN** the annotated table, plain ENTREZ list, and console summary are all produced
- **THEN** the annotated table contains one row per gene that passed the expression filter in at least one dataset
