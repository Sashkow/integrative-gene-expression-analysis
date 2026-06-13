## ADDED Requirements

### Requirement: Hard block-mask type in validate_imputation

`validate_imputation()` SHALL accept `mask_type = "gene_dataset_block_hard"` as a third masking strategy. This strategy SHALL build the (gene, dataset) pair pool from all genes present in a given dataset, without restricting to genes observed in all datasets. Any gene observed in ≥1 dataset is eligible for having one of its observed blocks masked.

#### Scenario: Hard masking includes partially-observed genes
- **WHEN** `validate_imputation()` is called with `mask_type = "gene_dataset_block_hard"` on a matrix where gene G is present in 3/6 datasets
- **THEN** gene G SHALL be included in the maskable pair pool for each of its 3 datasets

#### Scenario: Hard masking still respects min_obs_per_gene
- **WHEN** masking a block would leave a gene row with fewer than `min_obs_per_gene` observed cells
- **THEN** that block SHALL be skipped and the skip count SHALL be reported

### Requirement: Failure logging in validation results

`validate_imputation()` SHALL include two additional columns in its output data frame: `converged` (logical, TRUE when imputation succeeded, FALSE when tryCatch caught an error) and `error_message` (character, the error text on failure or NA on success). These columns SHALL be present for all mask types, not just hard masking.

#### Scenario: softImpute converges
- **WHEN** softImpute completes without error for a given repeat
- **THEN** the result row SHALL have `converged = TRUE` and `error_message = NA`

#### Scenario: softImpute fails (biScale or ALS divergence)
- **WHEN** softImpute throws an error during a repeat
- **THEN** the result row SHALL have `converged = FALSE`, accuracy metrics as `NA`, and `error_message` containing the error text from `conditionMessage(e)`

#### Scenario: Failure does not abort the validation loop
- **WHEN** softImpute fails on repeat 1 of 3
- **THEN** the function SHALL continue to repeats 2 and 3
