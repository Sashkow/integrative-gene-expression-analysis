## Why

The block-mask cross-validation (Table 3) restricts its gene pool to genes present in all 6 datasets. Masking one dataset's block from such a gene still leaves 5/6 observed — the easiest imputation target. Genes present in only 2–3 datasets are the ones that actually need imputation the most, yet their accuracy is never measured. The article already acknowledges this limitation (§Limitations). Running the existing pipeline with a hard masking strategy and varying `min_datasets` thresholds would either confirm that softImpute degrades gracefully or reveal a cliff.

## What Changes

- Add a new mask type `"gene_dataset_block_hard"` to the existing `validate_imputation()` function in `imputation.R` — same as `gene_dataset_block` but drops the restriction that genes must be present in all datasets; any gene with ≥2 dataset-observations is eligible for block masking
- Ensure softImpute failures in validation are logged (not fatal) — the `tryCatch` exists but verify it captures biScale/convergence failures cleanly and records failure counts
- Create 4 config YAML variants based on the 6ds article config, each with a different `min_datasets` threshold (2, 3, 4, 5) and `mask_type: "gene_dataset_block_hard"`
- Run via the existing `run_phase2b.R` pipeline with `--validate_only` flag, report per-threshold: Pearson r, RMSE, MAE, number of softImpute attempts, convergence failures, failure rate (%)

## Capabilities

### New Capabilities
- `hard-block-mask-type`: New `mask_type = "gene_dataset_block_hard"` in `validate_imputation()` that allows masking any gene regardless of coverage, plus failure-rate reporting in the validation output
- `hard-mask-configs`: 4 config variants for the 6ds dataset at `min_datasets` = 2, 3, 4, 5 with hard block masking enabled

### Modified Capabilities

## Impact

- `scripts/integrative_analysis/phase2b_direct_merge/imputation.R` — new mask type branch in `validate_imputation()`, possible improvements to tryCatch error logging
- New config YAMLs under `config_yehor_sashko/`
- Output to new subdirectories under `output/yehor_sashko/`
- No changes to the main pipeline logic, DE analysis, or existing results
