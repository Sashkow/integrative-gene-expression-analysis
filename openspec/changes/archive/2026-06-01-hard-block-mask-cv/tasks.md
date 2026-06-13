## 1. Hard block-mask type in validate_imputation

- [x] 1.1 Add `"gene_dataset_block_hard"` to `match.arg()` in `validate_imputation()` and implement the pair-pool construction — same as `gene_dataset_block` but without the `n_ds_per_gene == n_datasets_total` filter
- [x] 1.2 Add `converged` (logical) and `error_message` (character) columns to the results data frame for all mask types — populate from existing tryCatch, store `conditionMessage(e)` on failure
- [x] 1.3 Verify the existing tryCatch at line ~565 catches biScale and ALS convergence errors without aborting the repeat loop

## 2. Config variants

- [x] 2.1 Create 4 config YAMLs under `config_yehor_sashko/` based on the 6ds article config, each with `mask_type: "gene_dataset_block_hard"`, `min_obs_per_gene: 1`, and `min_datasets` set to 2, 3, 4, 5 respectively, with distinct output directories
- [x] 2.2 Verify each config is runnable with `Rscript run_phase2b.R --config <path> --validate_only`

## 3. Run and report (hard block mask)

- [x] 3.1 Execute all 4 configs with `--validate_only`, collect `imputation_validation.csv` from each output directory
- [x] 3.2 Inspect results: per-tier accuracy metrics and failure rates (converged vs total attempts)

## 4. Progressive-tax block masking

- [x] 4.1 Add `"progressive_tax_block"` mask type to `validate_imputation()` — for each repeat, compute `n_min` = minimum number of datasets any gene appears in across the matrix; then for every gene with coverage > `n_min`, mask all but `n_min` randomly selected datasets (genes already at `n_min` are left untouched). This equalizes all genes to worst-case coverage.
- [x] 4.2 Create config YAML(s) under `config_yehor_sashko/` with `mask_type: "progressive_tax_block"` for the 6ds dataset at `min_datasets: 2` (the lowest threshold that includes all coverage tiers)
- [x] 4.3 Run with `--validate_only`, inspect results — accuracy now reflects "what if every gene had the sparsity of the poorest gene"

## 5. Normalize masking fraction across all mask types

- [x] 5.1 Extract masking logic into separate functions (`mask_random_cells()`, `mask_gene_dataset_block()`, `mask_progressive_tax_block()`) that each return `list(mask_idx, n_blocks)` — keeps `validate_imputation()` clean
- [x] 5.2 Make all mask types respect `leave_out_fraction` so they each mask ~the same percentage of observed values. For progressive tax: instead of masking every gene down to n_min, randomly sample a subset of genes to mask down (or mask fewer blocks per gene) until the total masked cells ≈ `leave_out_fraction * n_observed`. For block-mask types: already controlled by `n_mask_target`, just verify.
- [x] 5.3 Rerun all configs with normalized masking fraction, compare results on equal footing

## 6. Hard block-mask with min_datasets=1

- [x] 6.1 Create config YAML with `min_datasets: 1`, `mask_type: "gene_dataset_block_hard"`, `min_obs_per_gene: 1` — includes all 17,531 genes; masking pool excludes genes present in only 1 dataset (they can't lose a block and retain any data), so the pool covers genes in 2–6 datasets while 1/6 genes sit in the matrix providing structural support
- [x] 6.2 Run with `--validate_only`, inspect results and failure rate
