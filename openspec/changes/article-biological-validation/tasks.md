## 1. Prior-study concordance analysis (re-run for 6ds)

- [x] 1.1 Write R script `scripts/integrative_analysis/article_validation/prior_study_concordance.R` — EXISTS, needs base_dir update from 7ds to 6ds
- [x] 1.2 Archive 7ds outputs: move `output/article_validation/prior_study_concordance/` contents to `output/article_validation/prior_study_concordance/archive/7ds/` (preserve clustered arrows and transition matrix outputs)
- [x] 1.3 Update `base_dir` in `prior_study_concordance.R` and all `volcano_arrows*.R` scripts from `phase2b_1_2_yehor_7ds_no_37653_enriched_sashko` to `phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko`
- [x] 1.4 Re-run `prior_study_concordance.R` to regenerate concordance table, scatter plot, volcano plot, and Venn diagram with 6ds DEGs (447 softImpute, 221 intersection, 231 gained)
- [x] 1.5 Re-run `volcano_arrows_clustered.R` and `transition_matrix.R` to regenerate clustered arrows and transition matrix/heatmap with 6ds data
- [x] 1.6 Verify updated concordance numbers: 231 gained, 211/231 matched (91.3%), 188/211 same direction (89.1%), 141 FDR<0.05 (66.8%), Pearson r=0.606

## 1b. Internal concordance: 6ds softImpute vs 6ds intersection-only

- [x] 1b.1 Write R script `internal_concordance.R` comparing softImpute vs intersection-only on 8,260 shared genes — scatter, volcano, clustered arrows, transition matrix
- [x] 1b.2 logFC scatter plot: r=0.9982, 97.8% same direction across 8,260 shared genes
- [x] 1b.3 Volcano plot of intersection-only full limma with softImpute-only DEGs highlighted (23 genes)
- [x] 1b.4 Clustered arrows plot: 8 clusters from 216 shared DEGs
- [x] 1b.5 Transition matrix/heatmap: 216 DEGs in both, 23 softImpute-only, 5 intersection-only

## 2. GO/KEGG enrichment analysis (re-run for 6ds)

- [x] 2.1 Write R script `scripts/integrative_analysis/article_validation/deg_enrichment.R` — EXISTS, needs base_dir update
- [ ] 2.2 Update `base_dir` in `deg_enrichment.R` from 7ds to 6ds path and update label strings (380→447, 182→221)
- [ ] 2.3 Re-run `deg_enrichment.R` to regenerate enrichment tables and dot plots with 6ds DEGs
- [ ] 2.4 Verify updated enrichment term counts and document

## 3. Mean-imputation benchmark

- [ ] 3.1 Add per-gene-mean and per-batch-mean imputation methods to the cross-validation framework in `scripts/integrative_analysis/phase2b_direct_merge/imputation.R`
- [ ] 3.2 Run cross-validation with mean-imputation baselines under both random-cell and block masking
- [ ] 3.3 Produce updated comparison table (Pearson r, RMSE, MAE) for all methods

## 4. ComBat sensitivity analysis

- [ ] 4.1 Write R script `scripts/integrative_analysis/article_validation/combat_sensitivity.R` that extracts ComBat batch parameters (gamma, delta) from intersection-only (8,260 genes) and full imputed (17,531 genes) runs using 6ds output
- [ ] 4.2 Compute and plot parameter correlation + report max deviation

## 5. Article text integration

- [ ] 5.1 Add Results subsection "Concordance with prior single-platform analysis" with prior-study concordance findings (using 6ds numbers), scatter plot reference, and volcano plot reference
- [ ] 5.2 Add Results subsection "Functional enrichment of differentially expressed genes" with enrichment findings and dot plot reference
- [ ] 5.3 Add Discussion paragraph interpreting biological validation results
- [ ] 5.4 Update imputation accuracy table with mean-imputation baseline rows
- [ ] 5.5 Add ComBat sensitivity results (as subsection or supplementary table reference)
- [ ] 5.6 Add Venn diagram figure and caption to Results
