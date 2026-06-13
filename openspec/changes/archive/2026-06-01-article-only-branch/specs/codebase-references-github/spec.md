## ADDED Requirements

### Requirement: GitHub permalink references in codebase_references.md
All file references in `codebase_references.md` SHALL use GitHub permalink format with the article branch commit SHA: `https://github.com/Sashkow/integrative-gene-expression-analysis/blob/<sha>/path/to/file`. Data files that exist in `articles/imputation_article/data/` SHALL be referenced via that path. Both scripts and data files SHALL be linked.

#### Scenario: Script references use GitHub links
- **WHEN** reading codebase_references.md
- **THEN** each script reference SHALL be a clickable GitHub permalink
- **AND** the SHA SHALL match the article branch commit

#### Scenario: Data references prefer article data path
- **WHEN** a data file exists in both `output/` and `articles/imputation_article/data/`
- **THEN** the reference SHALL prefer the `articles/imputation_article/data/` path
