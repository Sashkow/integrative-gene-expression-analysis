## Why

The repo (`github.com/Sashkow/integrative-gene-expression-analysis`) contains the full project: multiple articles, MCP servers, personal notes, conference materials, etc. Potential readers/reviewers of the imputation article (BMC Bioinformatics) shouldn't have to navigate this clutter. A clean branch with only article-related files makes the codebase reviewable and the results reproducible.

Additionally, `codebase_references.md` currently uses local paths. With an article-specific branch on GitHub, we can use permanent GitHub permalink references (commit SHA-based) to link article claims to the exact code and data that produced them.

## What Changes

- Create an orphan branch `article/imputation` containing only files related to the imputation article
- Install Git LFS and use it for large expression matrices (~150MB of `exprs_*.tsv`)
- Exclude unrelated directories (other articles, MCP servers, emotional_management, etc.)
- Exclude reference PDFs (73MB — readers can find via DOI)
- Include: article source (main.tex, figures, supplements, data/), pipeline scripts, configs, pipeline output (DE tables, validation, plots), article validation scripts and output
- Update `codebase_references.md` with GitHub permalink references to both scripts and data files on the branch
- Push the branch to origin

## Capabilities

### New Capabilities
- `article-branch-setup`: Create and populate the article-only orphan branch with Git LFS for large files
- `codebase-references-github`: Update codebase_references.md with GitHub permalink-style references

### Modified Capabilities

## Impact

- New branch `article/imputation` on GitHub
- Git LFS installed and configured for `*.tsv` files on that branch
- `articles/imputation_article/codebase_references.md` updated with GitHub links
- No changes to master branch
