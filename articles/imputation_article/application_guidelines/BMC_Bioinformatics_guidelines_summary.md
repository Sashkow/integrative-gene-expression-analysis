# BMC Bioinformatics -- Submission Guidelines Summary

Source: https://bmcbioinformatics.biomedcentral.com/submission-guidelines

## Journal Scope

BMC Bioinformatics is an open-access, peer-reviewed journal that considers articles
on computational algorithms, software, models and tools (including statistical methods,
ML, AI) for modelling and analysis of biological data, and systems biology.

**Article types:** Research Article, Software Article, Methodology Article, Database Article

This manuscript fits: **Research Article** (or possibly Methodology Article).

## Article Processing Charge

GBP 2,290 / USD 3,090 / EUR 2,590

## Required Sections (Research Article, in order)

### Front matter
1. **Title** -- concise, informative
2. **Authors** -- First Name, Middle Initial(s), Last Name; affiliations; corresponding author email
3. **Abstract** -- structured, max 350 words, sections: Background, Results, Conclusions. No references, minimal abbreviations.
4. **Keywords** -- 3-10 keywords

### Main body
5. **Background** -- context, aims, existing literature, why this study was necessary
6. **Results** -- findings including statistical results, tables, figures
7. **Discussion** -- implications in context of existing research, limitations
8. **Conclusions** -- main conclusions, importance and relevance
9. **Methods** -- clear description of all processes, interventions, comparisons; type of statistical analysis

Note: Results and Discussion can be combined. Methods typically comes last in BMC format.

### Back matter
10. **List of abbreviations** (if used)
11. **Declarations**
    - Ethics approval and consent to participate
    - Consent for publication
    - Availability of data and materials (MANDATORY -- all datasets must be available)
    - Competing interests
    - Funding
    - Authors' contributions
12. **Acknowledgements**
13. **References** -- Vancouver (numbered) style, via bmc-mathphys.bst
14. **Figure legends** (figures submitted separately)
15. **Tables**
16. **Additional files** (supplementary material, each with title and description)

## Formatting

- **LaTeX template:** `bmcart` document class (bmc_article_template.tex in this directory)
  - Alternative: Springer Nature `sn-jnl` class (sn-article_template.tex)
- **Figures:** submitted as separate files, not embedded in TeX. EPS, PDF, PNG, TIFF accepted. Min 300 DPI for raster.
- **Tables:** in the TeX file body, not as images.
- **References:** numbered Vancouver style. Use bmc-mathphys.bst or natbib with numbers.
- **Supplementary files:** each one described in "Additional Files" section with title and format description.

## Data and Code Availability

- **Strongly encouraged** to deposit datasets in publicly available repositories
- Code should be available (GitHub + Zenodo DOI recommended)
- Must include a "Availability of data and materials" declaration

## Key Differences from Current STAR Protocols Format

| Aspect | Current (STAR Protocols) | Required (BMC Bioinformatics) |
|--------|-------------------------|-------------------------------|
| Document class | article | bmcart or sn-jnl |
| Abstract | unstructured | structured: Background / Results / Conclusions, max 350 words |
| Section order | Key resources, Step-by-step, Expected outcomes | Background, Results, Discussion, Conclusions, Methods |
| Methods placement | middle of paper (Steps) | end of paper (before backmatter) |
| Figures | embedded with \includegraphics | legends only in TeX, files submitted separately |
| References | \begin{thebibliography} | .bib file + bmc-mathphys.bst (numbered, Vancouver) |
| Declarations | Resource availability only | Full set: ethics, competing interests, funding, contributions, data availability |
| Supplements | separate .tex documents | "Additional files" described in backmatter |

## Template Files in This Directory

- `bmc_article_template.tex` -- BMC article LaTeX template (v1.06)
- `bmcart.cls` -- BMC article document class
- `bmc_article.cls` -- alternative BMC class
- `bmc-mathphys.bst` -- BMC bibliography style
- `sn-jnl.cls` -- Springer Nature journal class (newer alternative)
- `sn-article_template.tex` -- Springer Nature article template
- `sn-mathphys.bst` -- Springer Nature bibliography style
