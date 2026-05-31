## ADDED Requirements

### Requirement: Orphan branch with article-only files
The branch `article/imputation` SHALL be an orphan branch containing only files related to the imputation article. It SHALL NOT contain unrelated project files (other articles, MCP servers, personal notes, conference materials).

#### Scenario: Branch contains article manuscript
- **WHEN** checking out `article/imputation`
- **THEN** `articles/imputation_article/main.tex` and `main.pdf` SHALL exist
- **AND** `articles/imputation_article/supplement1/` and `supplement2/` SHALL exist
- **AND** `articles/imputation_article/figures/` SHALL contain all referenced figures

#### Scenario: Branch contains pipeline scripts
- **WHEN** checking out `article/imputation`
- **THEN** `scripts/integrative_analysis/phase2b_direct_merge/` SHALL contain the pipeline scripts and configs
- **AND** `scripts/integrative_analysis/phase5_validation/` SHALL contain validation scripts and `config_validation_6ds.yaml`
- **AND** `scripts/integrative_analysis/article_validation/` SHALL contain concordance and enrichment scripts

#### Scenario: Branch contains pipeline output
- **WHEN** checking out `article/imputation`
- **THEN** DE tables, summary files, validation CSVs, and plots from the 6ds run SHALL be present
- **AND** expression matrices SHALL be tracked via Git LFS

#### Scenario: Branch excludes unrelated files
- **WHEN** checking out `article/imputation`
- **THEN** no MCP server directories, other articles, personal notes, or conference materials SHALL exist
- **AND** `articles/imputation_article/references/` (PDFs) SHALL NOT be included
