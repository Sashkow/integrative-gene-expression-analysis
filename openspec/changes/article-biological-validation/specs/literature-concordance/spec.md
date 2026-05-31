## ADDED Requirements

### Requirement: Literature concordance with published placental expression data

The validation MUST include a systematic PubMed search for published wet-lab data on the 447 DEGs. For each gene with published placental qPCR/RT-PCR/IHC/Western blot data reporting direction of change across trimesters, the concordance with our pipeline logFC MUST be recorded.

#### Scenario: Literature mining and curation

Given the 447 DEG gene symbols
When PubMed is searched per gene for placental expression studies (qPCR, RT-PCR, IHC, Western, protein)
And Human Protein Atlas placenta tissue data is checked
Then a curated CSV is produced with columns: gene, PMID, method, tissue, published_direction, our_logFC, concordant
And priority is given to well-known placental markers (LEP, FLT1, ENG, PAPPA, ADAM12, HSD3B1, CGA, CGB, GH2, PSG family) and to the ~226 gained DEGs

#### Scenario: Concordance summary for article

Given the curated literature concordance CSV
When concordance rate is computed (same-direction / total with published data)
Then a summary table is formatted for main.tex as a new Results subsection
And the Discussion interprets the concordance rate in context of pipeline validation
And special attention is given to gained DEGs with literature support (directly addresses the reviewer artefact concern)
