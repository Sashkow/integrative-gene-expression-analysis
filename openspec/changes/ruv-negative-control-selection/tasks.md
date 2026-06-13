## 1. Data preparation

- [ ] 1.1 Create script skeleton at `scripts/integrative_analysis/article_validation/select_ruv_empirical_controls.R` with library imports (limma, AnnotationDbi, org.Hs.eg.db, openxlsx) and path constants
- [ ] 1.2 Load GSE9984 and GSE100051 expression matrices from `data/mapped/yehor/8_yehor_preprocessed_datasets_enriched_sashko/`
- [ ] 1.3 Load sample metadata and split samples by trimester (1st vs 2nd) for each dataset

## 2. Expression filtering and variance ranking

- [ ] 2.1 Compute per-gene average expression across all samples within each dataset
- [ ] 2.2 Apply 25th-90th percentile expression filter per dataset
- [ ] 2.3 Compute CV (SD/mean) for each gene within each dataset
- [ ] 2.4 Print filter step counts (total genes → after expression filter → shared across both datasets)

## 3. DE exclusion

- [ ] 3.1 Run per-dataset limma DE (1T vs 2T) within GSE9984 and GSE100051 separately, exclude genes with FDR < 0.05
- [ ] 3.2 Load full pipeline DE results from `output/yehor_sashko/phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko/difexp_softimpute_combat_ref.tsv`, exclude genes with adj.P.Val < 0.05
- [ ] 3.3 Print counts after DE exclusion

## 4. Cross-reference annotations

- [ ] 4.1 Load Eisenberg-Levanon housekeeping list from `data/reference/integration_methods_references/ruv_housekeeping_genes_search/eisenberg_levanon_HK_genes.txt`, map symbols to ENTREZ, annotate candidates with housekeeping_flag
- [ ] 4.2 Curate literature DEG exclusion list: collect ENTREZ IDs of genes validated as 1_2 DEGs by qPCR/blot in placenta studies (Prater 2021 top validated genes, supplement 2 literature verdicts for unstable genes). Exclude from candidates.
- [ ] 4.3 Add gene symbols via org.Hs.eg.db mapping

## 5. Output

- [ ] 5.1 Build annotated per-gene table with all columns (ENTREZID, SYMBOL, mean expression, CV, ranks, DE status, housekeeping flag, literature DEG, final decision, exclusion reason)
- [ ] 5.2 Write annotated table to `output/article_validation/ruv_empirical_control_genes/ruv_empirical_candidates.xlsx`
- [ ] 5.3 Write plain ENTREZ ID list for final selected genes to `output/article_validation/ruv_empirical_control_genes/ruv_empirical_control_genes_1_2.txt`
- [ ] 5.4 Print console summary with gene counts at each step (funnel: total → expression-filtered → shared → non-DE → non-literature-DEG → final)

## 6. Verification

- [ ] 6.1 Run the script end-to-end and verify all output files are produced
- [ ] 6.2 Spot-check: verify ACTB and GAPDH are excluded (known unstable in placenta)
- [ ] 6.3 Spot-check: verify known stable housekeeping genes (YWHAZ, TBP, SDHA) are present if they pass empirical filters
- [ ] 6.4 Compare final list size against existing `build_ruv_control_genes.R` output to assess overlap
