## ADDED Requirements

### Requirement: Config variants for hard block-mask validation

Four config YAML files SHALL be created under `config_yehor_sashko/`, each based on the existing `config_phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko.yaml` with the following overrides:
- `validation.mask_type`: `"gene_dataset_block_hard"`
- `coverage.min_datasets`: 2, 3, 4, or 5 respectively
- `validation.min_obs_per_gene`: 1
- `paths.output`: a distinct output directory per threshold (e.g. `output/yehor_sashko/hard_block_mask_cv/min_ds_2/`)

Each config SHALL use the same 6 datasets, phenodata, imputation parameters, and normalization settings as the base config.

#### Scenario: Config for min_datasets=2
- **WHEN** the config with `min_datasets: 2` is loaded
- **THEN** the pipeline SHALL include all genes present in ≥2 of the 6 datasets and validate with hard block masking

#### Scenario: Config for min_datasets=5
- **WHEN** the config with `min_datasets: 5` is loaded
- **THEN** the pipeline SHALL include all genes present in ≥5 of the 6 datasets and validate with hard block masking

### Requirement: Runnable via existing pipeline

Each config SHALL be runnable via `Rscript run_phase2b.R --config <config_path> --validate_only` with no modifications to `run_phase2b.R` beyond what the hard-block-mask-type spec requires.

#### Scenario: Validate-only run with hard masking config
- **WHEN** `run_phase2b.R --validate_only` is invoked with one of the hard-mask configs
- **THEN** it SHALL run imputation validation, write `imputation_validation.csv` to the config's output directory, and exit without running DE analysis
