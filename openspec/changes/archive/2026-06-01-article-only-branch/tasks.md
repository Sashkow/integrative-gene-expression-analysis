## 1. Install Git LFS

- [x] 1.1 Install git-lfs (`sudo apt-get install git-lfs`) and run `git lfs install`

## 2. Create orphan branch and populate

- [x] 2.1 Create orphan branch `article/imputation` (`git checkout --orphan article/imputation`)
- [x] 2.2 Remove all auto-staged files (`git rm -rf --cached .`)
- [x] 2.3 Create article-branch `.gitignore` (allow articles/, output/, article_validation scripts; exclude references/, deconvolution data, build artifacts)
- [x] 2.4 Set up LFS tracking for large expression matrices
- [x] 2.5 Add pipeline scripts: `scripts/integrative_analysis/phase2b_direct_merge/`, `scripts/integrative_analysis/phase5_validation/`
- [x] 2.6 Add article validation scripts: `scripts/integrative_analysis/article_validation/`
- [x] 2.7 Add article manuscript: `articles/imputation_article/` (exclude `references/`, `data/deconvolution data/`, build artifacts)
- [x] 2.8 Add pipeline output: `output/yehor_sashko/phase2b_1_2_yehor_6ds_no_37653_no_22490_enriched_sashko/` (6 exprs_*.tsv via LFS)
- [x] 2.9 Add validation output: `output/yehor_sashko/phase5_validation_6ds/`
- [x] 2.10 Add article validation output: `output/article_validation/`
- [x] 2.11 Add balanced reference DE table: `output/yehor_sashko/phase2b_1_2_yehor_2ds_balanced/difexp_softimpute_combat.tsv`
- [x] 2.12 Add phenodata: `data/phenodata/samples.csv`, `data/phenodata/samples_evidence.md`
- [x] 2.13 Add `install_packages.R`, `CLAUDE.md`
- [x] 2.14 Commit the article branch (SHA: 5cc8489)

## 3. Update codebase_references.md with GitHub links

- [x] 3.1 Get the commit SHA from the article branch
- [x] 3.2 Rewrite `codebase_references.md` with GitHub permalink references (prefer `articles/imputation_article/data/` paths for data)
- [x] 3.3 Amend the commit with the updated codebase_references.md

## 4. Push

- [x] 4.1 Push `article/imputation` branch to origin (6 LFS objects, 153 MB uploaded)
- [x] 4.2 Switch back to master
