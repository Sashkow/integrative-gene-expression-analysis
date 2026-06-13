## Context

The current `validate_imputation()` in `imputation.R` supports `mask_type = "gene_dataset_block"` but restricts the gene pool to genes present in all N datasets (line 488: `n_ds_per_gene == n_datasets_total`). This is necessary because masking a block from a partially-observed gene can leave the row too sparse for `biScale()`, crashing the entire softImpute run.

The existing pipeline already has:
- `tryCatch` around imputation calls (line 565) that records `NA` metrics on failure
- `--validate_only` flag in `run_phase2b.R` to skip DE analysis
- `min_obs_per_gene` parameter to skip blocks that would leave a row too sparse
- Config-driven `mask_type` and `min_datasets` parameters

## Goals / Non-Goals

**Goals:**
- Add a `"gene_dataset_block_hard"` mask type that allows masking any gene with ≥2 dataset-observations
- Ensure softImpute failures are logged with error messages and counted, not swallowed silently
- Create 4 config variants at `min_datasets` 2–5 to run through the existing pipeline
- Report failure rates alongside accuracy metrics

**Non-Goals:**
- New runner script or validation function — reuse existing `validate_imputation()` and `run_phase2b.R`
- Modifying the production imputation path or existing DE results
- Fixing softImpute to handle sparse rows

## Decisions

**1. New mask type alongside existing ones**

Add `"gene_dataset_block_hard"` as a third option in `match.arg(mask_type, ...)`. The implementation is a copy of the `gene_dataset_block` branch with the `n_ds_per_gene == n_datasets_total` filter removed — eligible pairs are any (gene, dataset) where the gene is observed in that dataset, regardless of how many other datasets it appears in.

*Alternative considered*: Adding a flag like `restrict_to_full_coverage = FALSE` to the existing `gene_dataset_block` path. Rejected because it changes the semantics of an existing option — safer to keep the existing behavior untouched and add a new clearly-named variant.

**2. Failure logging in tryCatch**

The existing `tryCatch` at line 565 already catches errors and records `NA` metrics. Verify it also captures the error message in the output data frame (current code has `conditionMessage(e)` in cat output but doesn't store it in the results). Add an `error_message` column and a `converged` logical column to the results data frame.

**3. Config variants, not a sweep script**

Create 4 separate config YAMLs based on `config_phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko.yaml`, each differing only in:
- `coverage.min_datasets`: 2, 3, 4, or 5
- `validation.mask_type`: `"gene_dataset_block_hard"`
- `paths.output`: separate output directory per threshold

Run each with: `Rscript run_phase2b.R --config <config> --validate_only`

**4. min_obs_per_gene stays at default (or set to 1)**

For the stress test, set `min_obs_per_gene: 1` in the configs to allow maximum masking. Blocks that would leave a gene with 0 observations are still skipped (unavoidable).

## Risks / Trade-offs

**[softImpute crashes for low-coverage tiers]** → Expected and informative. The tryCatch ensures the pipeline doesn't abort; the failure is recorded as a data point.

**[biScale fails globally, not per-gene]** → A single ultra-sparse row can crash the whole run. With `min_obs_per_gene: 1`, rows with just 1 cell remaining after masking may cause this. If tier 2 fails 100% of the time, that's the finding.

**[Noisy metrics for small gene pools]** → Tier 5/6 may have fewer genes than tier 6/6. The output includes `n_genes_in_tier` / `n_masked_cells` so the reader can judge.

**[Runtime]** → 4 configs × `--validate_only` ≈ 4 × ~2-5 min = 8-20 min total. Acceptable.
