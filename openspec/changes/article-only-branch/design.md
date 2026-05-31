## Context

Repo is at `https://github.com/Sashkow/integrative-gene-expression-analysis` with a single `master` branch. The `.gitignore` excludes `articles/`, `output/`, `one_off_scripts/`, and large data dirs. The article-related files total ~330MB (of which ~150MB are expression matrices suitable for LFS, 73MB are reference PDFs to exclude).

## Goals / Non-Goals

**Goals:**
- Clean branch with only article-relevant files
- Large files tracked via Git LFS
- GitHub permalink references in codebase_references.md

**Non-Goals:**
- Modifying master branch
- Changing the pipeline code itself
- Including reference PDFs (available via DOI)

## Decisions

### 1. Orphan branch approach
Use `git checkout --orphan article/imputation` to create a branch with no shared history with master. This keeps the branch clean — no unrelated commits in history. Selectively copy/add the needed files.

### 2. Git LFS for expression matrices
Track `output/**/*.tsv` files >1MB via Git LFS. The 6 expression matrices (~150MB total) are the main LFS targets. Small TSVs (DE tables, validation results) stay in regular git.

### 3. File inventory for the branch

**From master (already tracked):**
- `scripts/integrative_analysis/phase2b_direct_merge/` — pipeline scripts + all configs
- `scripts/integrative_analysis/phase5_validation/` — validation scripts + configs
- `data/phenodata/samples.csv`, `samples_evidence.md`
- `install_packages.R`

**Force-add (currently gitignored):**
- `articles/imputation_article/` — main.tex, supplements, figures, data/logfc_stability/, codebase_references.md, xlsx files
  - EXCLUDE: `references/` (73MB PDFs), `data/deconvolution data/`, build artifacts
- `scripts/integrative_analysis/article_validation/` — concordance, enrichment, logfc_stability scripts
- `output/yehor_sashko/phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko/` — all pipeline output (LFS for large TSVs)
- `output/yehor_sashko/phase5_validation_6ds/` — validation TSVs + plots
- `output/article_validation/` — concordance + enrichment results
- `output/yehor_sashko/phase2b_1_2_yehor_2ds_balanced/difexp_softimpute_combat.tsv` — balanced reference DE table

**EXCLUDE from branch:**
- All other `articles/` directories
- `mcp-demo/`, `mcp-vectordb/`, `igea-vectordb-from-demo/`, `paper-search-mcp/`
- `emotional_management/`, `yehor/`, `presentation/`, `docs/`, `docs_to_ingest/`
- `GSE93520_redownload/`, `openspec/`, `notebooks/`, `tests/`
- `one_off_scripts/` (figures already baked into article dir)
- `articles/imputation_article/references/` (73MB)

### 4. GitHub permalink format
Use `https://github.com/Sashkow/integrative-gene-expression-analysis/blob/<sha>/path` with the commit SHA from the article branch. For data files that live in `articles/imputation_article/data/`, prefer that path over the output path.

## Risks / Trade-offs

- [Risk] Git LFS requires `sudo apt-get install git-lfs` → Mitigation: install as first step, verify before proceeding
- [Risk] GitHub free tier LFS limit (1GB storage) → Mitigation: ~150MB of TSVs fits well within limit
- [Risk] Orphan branch means no shared history → Mitigation: intentional — keeps article branch reviewable without 100+ unrelated commits
