# Codebase References for Article Numbers

Maps article claims to the code and data that produced them.

GitHub base: `https://github.com/Sashkow/integrative-gene-expression-analysis`
Branch: `article/imputation`

[B]: https://github.com/Sashkow/integrative-gene-expression-analysis/blob/article/imputation

## Source data

- **6ds pipeline config**: [`config_phase2b_1_2_yehor_6ds_...yaml`][B]/scripts/integrative_analysis/phase2b_direct_merge/config_yehor_sashko/config_phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko.yaml)
- **6ds pipeline output**: [`output/yehor_sashko/phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko/`][B]/output/yehor_sashko/phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko)
- **6ds validation config**: [`config_validation_6ds.yaml`][B]/scripts/integrative_analysis/phase5_validation/config_validation_6ds.yaml)
- **6ds validation output**: [`output/yehor_sashko/phase5_validation_6ds/`][B]/output/yehor_sashko/phase5_validation_6ds)
- **Balanced 2ds reference**: [`difexp_softimpute_combat.tsv`][B]/output/yehor_sashko/phase2b_1_2_yehor_2ds_balanced/difexp_softimpute_combat.tsv)

## Results section

### Gene recovery (Table 1)
- Data: [`imputed_gene_group_coverage.csv`][B]/output/yehor_sashko/phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko/imputed_gene_group_coverage.csv) + [`summary.txt`][B]/output/yehor_sashko/phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko/summary.txt)
- Key numbers: 8,260 intersection, 17,531 at k=1, 29.4% missing

### Imputation accuracy (Table 3)
- Data: [`imputation_validation.csv`][B]/output/yehor_sashko/phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko/imputation_validation.csv)
- Block-mask: r=0.816, RMSE=1.523, MAE=1.165 (mean of 3 reps)
- Random-cell: from 7ds run (not yet re-run for 6ds) — r=0.994, RMSE=0.330, MAE=0.220

### DE counts (Table 4)
- Data: [`method_comparison.csv`][B]/output/yehor_sashko/phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko/method_comparison.csv)
- softImpute+ComBat-ref: 447 DEGs (369 up, 78 down)
- none+ComBat-ref: 221 DEGs (174 up, 47 down)
- softImpute+ComBat: 416 (346 up, 70 down)
- none+ComBat: 221 (175 up, 46 down)

### logFC stability (softImpute vs intersection-only, shared 8,260 genes)
- Script: [`logfc_stability.R`][B]/scripts/integrative_analysis/article_validation/logfc_stability.R)
- Data: [`logfc_stability.csv`][B]/articles/imputation_article/data/logfc_stability/logfc_stability.csv), [`logfc_stability.txt`][B]/articles/imputation_article/data/logfc_stability/logfc_stability.txt)
- r=0.9982, median |diff|=0.014, P99=0.086

### Subsampling validation
- Script: [`test1_first_trim_subsample.R`][B]/scripts/integrative_analysis/phase5_validation/test1_first_trim_subsample.R)
- Data: [`test1_first_trim_subsample.tsv`][B]/output/yehor_sashko/phase5_validation_6ds/test1_first_trim_subsample.tsv)
- CCC ≥ 0.88 at N=10, Jaccard 0.54 at N=10
- Retention 90% at N=10, 93% at N=30

### Balanced reference comparison
- Script: [`test1b_vs_balanced.R`][B]/scripts/integrative_analysis/phase5_validation/test1b_vs_balanced.R)
- Data: [`test1b_vs_balanced.tsv`][B]/output/yehor_sashko/phase5_validation_6ds/test1b_vs_balanced.tsv) vs [`difexp_softimpute_combat.tsv`][B]/output/yehor_sashko/phase2b_1_2_yehor_2ds_balanced/difexp_softimpute_combat.tsv) (416 DEGs, filtered at |logFC|>1, FDR<0.05)
- Shared: 409, Full-only: 38, Balanced-only: 7
- 91.5% (409/447) of full DEGs in balanced
- 98.3% (409/416) of balanced DEGs in full
- Jaccard: 0.901
- Note: test1b `overlap_vs_balanced` = retention of balanced DEGs in subsample (NOT % of full in balanced)

### Split-half
- Script: [`test3_split_half.R`][B]/scripts/integrative_analysis/phase5_validation/test3_split_half.R)
- Data: [`test3_split_half.tsv`][B]/output/yehor_sashko/phase5_validation_6ds/test3_split_half.tsv)
- Jaccard between halves: 0.40, CCC: 0.70
- Retention of full-run DEGs per half: 84%
- FDR-only retention: 51%, FDR-only DEGs: ~6,685

## Figures
- **PCA + staircase**: generated from pipeline output, final versions in [`figures/`][B]/articles/imputation_article/figures)
- **Validation plots**: [`plot_validation_results.R`][B]/scripts/integrative_analysis/phase5_validation/plot_validation_results.R) with [`config_validation_6ds.yaml`][B]/scripts/integrative_analysis/phase5_validation/config_validation_6ds.yaml), output in [`plots/`][B]/output/yehor_sashko/phase5_validation_6ds/plots)

## Prior-study concordance (article_validation)
- Script: [`prior_study_concordance.R`][B]/scripts/integrative_analysis/article_validation/prior_study_concordance.R)
- Data: [`gained_deg_concordance.csv`][B]/output/article_validation/prior_study_concordance/gained_deg_concordance.csv), [`concordance_summary.txt`][B]/output/article_validation/prior_study_concordance/concordance_summary.txt)
- 231 gained DEGs, 89.1% same direction vs 2021, Pearson r=0.606

## Internal concordance (softImpute vs intersection-only)
- Script: [`internal_concordance.R`][B]/scripts/integrative_analysis/article_validation/internal_concordance.R)
- Data: [`transition_matrix_internal.csv`][B]/output/article_validation/internal_concordance/transition_matrix_internal.csv)
- r=0.9982 on 8,260 shared genes, 216 DEGs in both, 23 softImpute-only, 5 intersection-only

## Pipeline (main scripts)
- [`run_phase2b.R`][B]/scripts/integrative_analysis/phase2b_direct_merge/run_phase2b.R) — main pipeline runner
- [`imputation.R`][B]/scripts/integrative_analysis/phase2b_direct_merge/imputation.R) — imputation methods
- [`normalization.R`][B]/scripts/integrative_analysis/phase2b_direct_merge/normalization.R) — batch correction
- [`plot_na_staircase.R`][B]/scripts/integrative_analysis/phase2b_direct_merge/plot_na_staircase.R) — staircase visualization
