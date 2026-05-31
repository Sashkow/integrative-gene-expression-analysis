# Samples.csv Evidence Log

Source-of-truth for gestational age (GA), diagnosis, and tissue data in `data/phenodata/samples.csv`.
When modifying samples.csv, add a dated entry here documenting what changed and why.

## GA Evidence by Dataset

### GSE100051 — Soncin 2018
- **GA type**: 49 exact (weeks 4–16), 5 Term category-only
- **Source**: GEO sample characteristics had GA as text ("Placenta; gestational age: Week 8"). Parsed to integer weeks.
- **Term samples** (GSM2668750–54): GA originally "Term", cleared to empty. No exact GA in article or GEO. Falls back to trimester-based estimation (37–41).
- **Article**: `articles/references/dataset_articles/GSE100051_Soncin_2018.pdf`
- **Diagnosis**: All Healthy (article: normal placentas across gestation)
- **Modified**: 2026-05-10 — cleaned text GA to integers, cleared Term GA values

### GSE107824 — Soncin 2018
- **GA type**: 11 exact (weeks 8–39)
- **Source**: GEO sample titles encode GA (e.g., "CTB_8wks", "CTB_39wks")
- **GEO page**: `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107824`
- **Article**: `articles/references/dataset_articles/GSE107824_Soncin_2018.pdf`
- **Diagnosis**: All Healthy
- **Modified**: 2026-05-10 — set exact GA from GEO titles

### GSE122214 — Zhao 2019
- **GA type**: 4 with bounds (7–8 weeks)
- **Source**: Article Methods: "30–35 days after embryo transfer" = gestational weeks 7–8
- **Article**: `articles/references/dataset_articles/GSE122214_Zhao_2019.pdf`
- **Diagnosis**: All Healthy (4 spontaneously conceived controls)
- **Modified**: 2026-05-10 — set bounds 7–8 from article

### GSE18044 — Bruchova 2010
- **GA type**: 63 exact (weeks 35–42), 1 with bounds (35–42)
- **Source**: GEO series matrix per-sample characteristic "gestational age (weeks): N"
- **GEO data**: `data/raws/GSE18044_series_matrix.txt.gz`
- **Bounds source**: Article Table 1 — population GA median 39 (range 35–42). Sample P163 had "not available" on GEO.
- **Article**: `articles/references/dataset_articles/GSE18044_Bruchova_2010.pdf`
- **Diagnosis**: 64 Healthy (non-smokers), 12 Healthy, Smoking. Source: GEO series matrix titles ("non-smoker N", "heavy smoker N", "light smoker N")
- **Modified**: 2026-05-10 — set exact GA from series matrix, diagnosis from titles

### GSE22490 — Rull 2013
- **GA type**: 6 exact (weeks 4–13)
- **Source**: Original samples.csv (pre-existing data)
- **Article**: `articles/references/dataset_articles/GSE22490_Rull_2013.pdf`
- **Note**: Same research group as GSE37901 (Uuskula 2012). GSE22490 = first trimester (n=6), GSE37901 = second trimester (n=4).
- **Diagnosis**: All Healthy (original data)

### GSE27272 — Votavova 2011
- **GA type**: 37 exact for placenta subset (weeks 37–41), 183 total across all tissues
- **Source**: GEO series matrix per-sample characteristic "gestational age (weeks): N"
- **GEO data**: `data/raws/GSE27272_series_matrix.txt.gz`
- **Diagnosis**: 128 Healthy (non-smokers), 55 Healthy, Smoking. Source: GEO series matrix characteristic "smoking status: smoker/non-smoker"
- **Article reference**: Votavova et al., Placenta 2011
- **Modified**: 2026-05-10 — set exact GA from series matrix, diagnosis from smoking status

### GSE28551 — Sitras 2012
- **GA type**: 37 with bounds (T1: 9–12, Term: 38–40)
- **Source**: Article — T1: mean 71.2 ± 8 days = weeks 9–12; T3: mean 275 ± 8 days = weeks 38–40
- **Article**: `articles/references/dataset_articles/GSE28551_Sitras_2012.pdf`
- **Diagnosis**: All Healthy (original data)
- **Modified**: 2026-05-10 — set bounds from article

### GSE35574 — Guo 2013
- **GA type**: 40 exact (weeks 29–40, healthy subset)
- **Source**: Original samples.csv (pre-existing data)
- **Article**: `articles/references/dataset_articles/GSE35574_Guo_2013.pdf`
- **Diagnosis**: 40 Healthy, 54 not healthy (original data)

### GSE37901 — Uuskula 2012
- **GA type**: 4 exact (weeks 17, 17, 18, 19)
- **Source**: Article Methods — gestational days 120, 121, 126, 132 for the 4 second-trimester samples. GEO titles encoded decimal weeks (17.1, 17.3, 18, 18.8). Excel had corrupted decimals to Ukrainian date strings ("17.січ.", "17.бер.", "18.Сер").
- **Article**: `articles/references/dataset_articles/GSE37901_Uuskula_2012.pdf` (page 4: "second (n=4; gestational days 120, 121, 126, 132) trimester placentae")
- **GEO page**: `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37901` (titles: MPG_18, MPG_17.1, MPG_17.3, MPG_18.8)
- **Diagnosis**: All Healthy (original data)
- **Modified**: 2026-05-10 — fixed corrupted GA values using article gestational days

### GSE43942 — No article downloaded
- **GA type**: 7 healthy placenta samples, category-only (Term)
- **Source**: Original samples.csv trimester category. No exact GA or bounds found.
- **GEO series matrix**: Not downloaded. GEO page not checked for per-sample GA.
- **Diagnosis**: 7 Healthy, 5 not healthy (original data)
- **Status**: Needs GEO check for per-sample GA

### GSE47187 — No article downloaded
- **GA type**: 5 healthy placenta samples, category-only (Late Preterm)
- **Source**: Original samples.csv trimester category. No exact GA or bounds found.
- **GEO series matrix**: Not downloaded. GEO page not checked for per-sample GA.
- **Diagnosis**: 5 Healthy, 5 not healthy (original data)
- **Status**: Needs GEO check for per-sample GA

### GSE55439 — Wolfe 2014
- **GA type**: 24 with bounds (T1: 10–12, Term: 38–40)
- **Source**: Article — T1 collected at 10–12 weeks, Term at 38–40 weeks
- **Article**: `articles/references/dataset_articles/GSE55439_Wolfe_2014.pdf`
- **Note**: 24 arrays = 6 biological samples × 4 technical replicates (RNAlater preservation). Collapsed to 6 in plots via `technical_replicate_block`.
- **Diagnosis**: All Healthy (original data)
- **Modified**: 2026-05-10 — set bounds from article

### GSE6573 — Herse 2007
- **GA type**: 1 exact (week 39, healthy placenta subset)
- **Source**: Original samples.csv (pre-existing data)
- **Article**: `articles/references/dataset_articles/GSE6573_Herse_2007.pdf`
- **Diagnosis**: 2 Healthy, 2 not healthy (original data). 2 excluded.

### GSE73374 — Martin 2015
- **GA type**: 17 exact (weeks 35–41, healthy subset)
- **Source**: Original samples.csv (pre-existing data)
- **Article**: `articles/references/dataset_articles/GSE73374_Martin_2015.pdf`
- **Diagnosis**: 17 Healthy, 19 not healthy (original data)

### GSE73685 — Bukowski 2017
- **GA type**: 18 healthy placenta samples, category-only (Late Preterm, Term)
- **Source**: Original samples.csv trimester category. Article and GEO checked — no per-sample GA available, only term/preterm categories.
- **Article**: `articles/references/dataset_articles/GSE73685_Bukowski_2017.pdf`
- **GEO data**: `data/raws/GSE73685_series_matrix.txt.gz`
- **Diagnosis**: 144 Healthy, 31 not healthy, 8 excluded (original data)
- **Status**: No improvement possible — article only has term/preterm groupings

### GSE7434 — Huuskonen 2008
- **GA type**: 5 healthy placenta samples, category-only (Term)
- **Source**: Original samples.csv trimester category. Article says "full-term placentas" with "normal delivery and healthy newborns" — no specific GA weeks.
- **Article**: `articles/references/dataset_articles/GSE7434_Huuskonen_2008.pdf`
- **GEO data**: `data/raws/GSE7434_series_matrix.txt.gz` — characteristics say "Full term placenta" only
- **Diagnosis**: 5 Healthy (non-smokers), 5 Healthy, Smoking. Source: GEO series matrix titles ("JR-N placenta, smoking/non-smoking, XX/XY")
- **Modified**: 2026-05-10 — set diagnosis from GEO titles
- **Status**: No GA improvement possible

### GSE74341 — No article downloaded
- **GA type**: 21 healthy placenta samples, category-only (Late Preterm, Term)
- **Source**: Original samples.csv trimester category. No exact GA or bounds found.
- **GEO series matrix**: Not downloaded. GEO page not checked for per-sample GA.
- **Diagnosis**: 21 Healthy, 29 not healthy (original data)
- **Status**: Needs GEO check for per-sample GA

### GSE93520 — Zhang 2017
- **GA type**: 36 exact (weeks 6–10)
- **Source**: Original samples.csv (pre-existing data)
- **Article**: `articles/references/dataset_articles/GSE93520_Zhang_2017.pdf`
- **Diagnosis**: All Healthy (original data)

### GSE9984 — Mikheev 2008
- **GA type**: 12 with bounds (T1: 6–8, T2: 15–16, Term: 37–41)
- **Source**: Article — T1: 45–59 days (6–8 weeks), T2: 109–115 days (15–16 weeks), Term: C-section deliveries (standard term 37–41)
- **Article**: `articles/references/dataset_articles/GSE9984_Mikheev_2008.pdf`
- **Diagnosis**: All Healthy (original data)
- **Modified**: 2026-05-10 — set bounds from article

### GDS2080 / GSE4707 — Nishizawa 2007
- **GA type**: Not set in samples.csv. Article Table 1: Normal 32.9 ± 5.6 weeks (n=4 microarray), Severe PE 32.9 ± 4.0 weeks (n=10 microarray; early onset <31 wks, late onset ≥31 wks)
- **Source**: Article Table 1
- **Article**: `articles/references/dataset_articles/GSE4707_Nishizawa_2007.pdf`
- **Tissue**: "A central area of chorionic tissue was dissected, and the maternal deciduas and amnionic membranes were removed. We then dissected 1-cm-thick sections of placental vili from the central area between basal and chorionic plates."
- **Diagnosis**: 4 Normal, 10+ Severe PE (early and late onset). All cesarean, no labor.
- **Platform**: Agilent Whole Human Genome Oligo Microarray, two-color

### GSE10588 — Sitras 2009
- **GA type**: Not set in samples.csv. Article Table 1: PE 238 ± 25 days (~34 wks), Controls 277 ± 9 days (~39.6 wks)
- **Source**: Article Table 1
- **Article**: `articles/references/dataset_articles/GSE10588_Sitras_2009.pdf`
- **Tissue**: "Chorionic tissue was dissected from a standardised location — approximately 2 cm beside the umbilical cord insertion, from the middle layer of placenta midway between maternal and fetal surfaces." ~2 cm³ specimen, washed with saline.
- **Diagnosis**: 16 Severe PE, 21 Healthy controls. 11/16 PE delivered prematurely (<37 wks). 3 had HELLP.
- **Platform**: Applied Biosystems Human Genome Survey Microarray v.2.0

### GSE12216 — Sitras 2009
- **GA type**: Not set in samples.csv. Article Table 1: IUGR 236 ± 24 days (~33.7 wks), Controls 274 ± 5 days (~39.1 wks)
- **Source**: Article Table 1
- **Article**: `articles/references/dataset_articles/GSE12216_Sitras_2009.pdf`
- **Tissue**: Same protocol as GSE10588 — "Chorionic tissue was dissected from a standardized location (approximately 2 cm beside the umbilical cord insertion, from the middle layer of placenta midway between maternal and fetal surfaces)." ~2 cm³, washed.
- **Diagnosis**: 8 IUGR (4 also with PE), 8 Controls. GA assigned by ultrasound at 18–20 weeks.
- **Platform**: Applied Biosystems Human Genome Survey Microarray v.2.0

### GSE12767 — Founds 2009
- **GA type**: Not set in samples.csv. Article Table 1: CVS at 11.3 ± 0.6 weeks (controls), 11.4 ± 0.7 weeks (PE). Range 10.7–12.4 weeks.
- **Source**: Article Table 1
- **Article**: `articles/references/dataset_articles/GSE12767_Founds_2009.pdf`
- **Tissue**: Surplus CVS specimens — "villi grossly free of decidua and maternal blood were removed from the Amniomax in the Petri dish, placed in an Eppendorf tube, and snap-frozen in less than 10 min from CVS aspiration."
- **Diagnosis**: 4 PE (developed later), 8 Controls. Matched by parity, GA at CVS, and race.
- **Platform**: Affymetrix HG-U133 Plus 2.0
- **Modified**: 2026-05-10 — changed Biological.Specimen from "Chorion" to "Chorionic Villus Sampling"

### GSE13155 — Cox 2009
- **GA type**: Term (~38 weeks), 2 placentas. Not set in samples.csv.
- **Source**: Article: "Two normal placentas at term (~38 weeks)" from elective cesarean
- **Article**: `articles/references/dataset_articles/GSE13155_Cox_2009.pdf`
- **Tissue**: Human villous trees dissected from term placentas. Also mouse labyrinth (E17.5).
- **Diagnosis**: Normal only (no disease comparison). 8 human villi samples from 2 individuals.
- **Platform**: Affymetrix HG-U133 Plus 2.0

### GSE24129 — Nishizawa 2011
- **GA type**: Not set in samples.csv. Article Table 1: Control 38.1 ± 0.8 wks, FGR 37.3 ± 1.0 wks, PE 34.4 ± 1.8 wks
- **Source**: Article Table 1
- **Article**: `articles/references/dataset_articles/GSE24129_Nishizawa_2011.pdf`
- **Tissue**: "A central area of chorionic tissue was dissected, and the maternal deciduas and amnionic membranes were removed. We then dissected 1 cm sections of placental villi from four different central areas between the basal and chorionic plates."
- **Diagnosis**: 8 Severe PE, 8 FGR, 8 Controls. All cesarean, no labor.
- **Platform**: Affymetrix Human Exon 1.0 ST Array

### GSE30186 — Meng 2012
- **GA type**: Not set in samples.csv. Article: Control 273 ± 5 days (~39.0 wks), PE 255 ± 6 days (~36.4 wks)
- **Source**: Article
- **Article**: `articles/references/dataset_articles/GSE30186_Meng_2012.pdf`
- **Tissue**: "Tissue blocks (approximately 1 cm³ each) were dissected from the standard locations on the maternal face of the placentas... Villous portions were harvested by dissecting free of blood vessels and connective tissue and washing off adherent blood clots."
- **Diagnosis**: 6 PE, 6 Controls. Elective cesarean, no labor.
- **Platform**: Illumina HumanHT-12 V4 BeadChip

### GSE36083 — Pantham 2012
- **GA type**: 8–8.5 weeks (first trimester explant study). Not set in samples.csv.
- **Source**: Article: "Three first-trimester placenta (8–8.5 weeks gestation)"
- **Article**: `articles/references/dataset_articles/GSE36083_Pantham_2012.pdf`
- **Tissue**: First-trimester placentas "dissected into explants of approximately 400 mg wet weight and cultured with aPL ID2 (25 μg/mL) for 16 h." In vitro treatment study, not clinical.
- **Diagnosis**: In vitro — 3 treated with antiphospholipid antibody vs. 3 untreated. Placentas from elective termination.
- **Platform**: Affymetrix HG-U133 Plus 2.0

### GSE44711 — Blair 2013
- **GA type**: EOPET <34 weeks, controls GA-matched preterm. Majority 31–34 weeks. Not set in samples.csv.
- **Source**: Article: EOPET defined as diagnosis <34 weeks, controls GA-matched
- **Article**: `articles/references/dataset_articles/GSE44711_Blair_2013.pdf`
- **Tissue**: "Whole chorionic villi were sampled from placentas... Chorionic villi from at least two sites (centre and perimeter) on the fetal side of the placenta were sampled, rinsed of maternal blood."
- **Diagnosis**: 8 EOPET + 8 GA-matched normotensive preterm controls (expression subset from larger methylation cohort of 20+20).
- **Platform**: Illumina HT-12v4 Expression BeadChip

### GSE54618 — Jebbink 2015
- **GA type**: Not set in samples.csv. Cohort A (microarray): normotensive median 31+0 wks (range 27+0–38+4), PE median 32+3 wks (range 28+5–38+2). Microarray on 10+10 subset.
- **Source**: Article
- **Article**: `articles/references/dataset_articles/GSE54618_Jebbink_2015.pdf`
- **Tissue**: "Placental biopsies from a macroscopically viable (non-infarcted) central cotyledon from the maternal side were obtained immediately after delivery and stored in RNAlater." From PANDA biobank.
- **Diagnosis**: 10 PE + 10 normotensive (microarray subset from Cohort A of 14 PE + 17 normotensive).
- **Platform**: Illumina HumanHT-12 v4 Expression BeadChip

### GSE57050 — No article downloaded
- **GA type**: Not set in samples.csv
- **GEO source**: characteristics: "tissue: villous tissue"
- **Tissue**: Villous tissue (chorionic villi)
- **Status**: Needs article download for GA and diagnosis details

## Non-placenta datasets with GA modifications

### GSE70102 — Bianco 2016
- **GA type**: 20 with bounds (Healthy: 12–20, T13: 13–19, T18: 14–20, T21: 18–22)
- **Source**: Article Methods — per-group GA ranges
- **Article**: `articles/references/dataset_articles/GSE70102_Bianco_2016.pdf`
- **Modified**: 2026-05-10 — refined bounds from 13–22 to per-group ranges

### GSE91189 — Garrido-Gomez 2017
- **GA type**: 4 Late Preterm with bounds (25–34), 4 Term cat-only, 4 Second Trimester cat-only
- **Source**: Article — "5 sPE and 5 preterm labor cases with no signs of infection (nPTB) (25-34 wks)"
- **Article**: `articles/references/dataset_articles/GSE91189_Garrido-Gomez_2017.pdf`
- **GEO data**: `data/raws/GSE91189_series_matrix.txt.gz` — no GA in characteristics
- **Modified**: 2026-05-10 — set Late Preterm bounds 25–34 from article

## Tissue / Extraction Protocol Evidence

Key distinction: **Chorion** (chorion laeve) = fetal membrane layer, collected separately from placental disc. **Chorionic villi** = villous tree of the placenta. **Placenta** = generic label, often full-thickness biopsy or unspecified villous tissue.

### Datasets currently labeled "Chorion" — review of actual tissue

#### GSE12767 — Founds 2009
- **Current label**: Chorionic Villus Sampling
- **Actual tissue**: CVS biopsies (chorionic villus sampling specimens)
- **GEO source**: "snap frozen, banked CVS specimens"
- **Article**: Founds SA et al. "Altered global gene expression in first trimester placentas of women destined to develop preeclampsia." Placenta 2009;30:15–24
- **Modified**: 2026-05-10 — changed from "Chorion" to "Chorionic Villus Sampling"

#### GSE13155 — Cox 2009
- **Current label**: Chorionic villi (changed from Chorion)
- **Actual tissue**: Chorionic villi
- **GEO source**: source_name = "placental villus tree"; characteristics = "Normal healthy term placenta"
- **Article**: `articles/references/dataset_articles/GSE13155_Cox_2009.pdf` — "human villous trees dissected from term placentas"
- **Verdict**: Relabeled to "Chorionic villi"

#### GSE24129 — Nishizawa 2011
- **Current label**: Chorionic villi (changed from Chorion)
- **Actual tissue**: Chorionic villi
- **GEO source**: source_name = "Human placenta"; extraction protocol: "chorionic villous tissues"
- **Article**: `articles/references/dataset_articles/GSE24129_Nishizawa_2011.pdf` — "dissected 1 cm sections of placental villi from four different central areas between the basal and chorionic plates"
- **Verdict**: Relabeled to "Chorionic villi"

#### GSE36083 — Pantham 2012
- **Current label**: Chorionic villi (changed from Chorion)
- **Actual tissue**: Chorionic villi (first-trimester placental explants)
- **GEO source**: "First-trimester human placenta, untreated"; 8–8.5 weeks gestation, dissected into explants. At 8 weeks this is villous tissue.
- **Article**: `articles/references/dataset_articles/GSE36083_Pantham_2012.pdf` — "dissected into explants of approximately 400 mg wet weight"
- **Verdict**: Relabeled to "Chorionic villi"

#### GSE44711 — Blair 2013
- **Current label**: Chorionic villi (changed from Chorion)
- **Actual tissue**: Chorionic villi
- **GEO source**: source_name = "ChorionicVilli"; characteristics: "tissue: Chorionic Villi", "gestational age (weeks): 31.7"
- **Article**: `articles/references/dataset_articles/GSE44711_Blair_2013.pdf` — "Whole chorionic villi were sampled from placentas... from at least two sites (centre and perimeter) on the fetal side"
- **Verdict**: Relabeled to "Chorionic villi"

#### GSE57050
- **Current label**: Chorionic villi (changed from Chorion)
- **Actual tissue**: Chorionic villi (villous tissue)
- **GEO source**: characteristics: "tissue: villous tissue"
- **Article**: Not downloaded
- **Verdict**: Relabeled to "Chorionic villi"

#### GDS2080 / GSE4707 — Nishizawa 2007
- **Current label**: Chorionic villi (changed from Chorion)
- **Actual tissue**: Chorionic villi (confirmed by article)
- **GEO source**: "Placental tissue from normal pregnancy", Caesarean section biopsies, PE study
- **Article**: `articles/references/dataset_articles/GSE4707_Nishizawa_2007.pdf` — "dissected 1-cm-thick sections of placental vili from the central area between basal and chorionic plates"
- **Verdict**: Relabeled to "Chorionic villi"

#### GSE10588 — Sitras 2009
- **Current label**: Chorionic villi (changed from Chorion)
- **Actual tissue**: Chorionic villi (confirmed by article — same protocol as GSE12216)
- **GEO source**: "total RNA from placenta tissue"
- **Article**: `articles/references/dataset_articles/GSE10588_Sitras_2009.pdf` — "Chorionic tissue was dissected from a standardised location — approximately 2 cm beside the umbilical cord insertion, from the middle layer of placenta midway between maternal and fetal surfaces"
- **Verdict**: Relabeled to "Chorionic villi"

#### GSE12216 — Sitras 2009
- **Current label**: Chorionic villi (changed from Chorion)
- **Actual tissue**: Chorionic villi (confirmed by article — same protocol as GSE10588)
- **GEO source**: source_name = "placenta"
- **Article**: `articles/references/dataset_articles/GSE12216_Sitras_2009.pdf` — "Chorionic tissue was dissected from a standardized location (approximately 2 cm beside the umbilical cord insertion, from the middle layer of placenta midway between maternal and fetal surfaces)"
- **Verdict**: Relabeled to "Chorionic villi"

#### GSE30186 — Meng 2012
- **Current label**: Chorionic villi (changed from Chorion)
- **Actual tissue**: Chorionic villi (confirmed by article)
- **GEO source**: "placentas from preeclamptic pregnancies"
- **Article**: `articles/references/dataset_articles/GSE30186_Meng_2012.pdf` — "Villous portions were harvested by dissecting free of blood vessels and connective tissue and washing off adherent blood clots"
- **Verdict**: Relabeled to "Chorionic villi"

#### GSE54618 — Jebbink 2015
- **Current label**: Chorionic villi (changed from Chorion)
- **Actual tissue**: Placental cotyledon biopsies (maternal side)
- **GEO source**: "placental tissue biopsies"
- **Article**: `articles/references/dataset_articles/GSE54618_Jebbink_2015.pdf` — "Placental biopsies from a macroscopically viable (non-infarcted) central cotyledon from the maternal side were obtained immediately after delivery"
- **Verdict**: **Mislabeled** — should be "Placenta" or "Chorionic villi" (cotyledon biopsy from maternal side)

#### GSE73685 — Bukowski 2017
- **Current label**: Chorion (24 samples) + Placenta (21 samples) + other tissues
- **Actual tissue**: Chorion membrane (chorion laeve), collected separately from placental disc
- **Article**: `articles/references/dataset_articles/GSE73685_Bukowski_2017.pdf`
- **Verdict**: **Correctly labeled** — this is the only dataset where "Chorion" is accurate

### Datasets labeled "Placenta" — extraction protocol notes

#### GSE100051 — Soncin 2018
- **Tissue**: Whole placental tissue (villous tissue, mixed cell population)
- **Article quote**: "we investigated gene expression in a random sampling of villous tissue, containing a mixed cell population, with stromal tissues as well as trophoblasts"
- **Verdict**: "Placenta" acceptable — whole villous tissue

#### GSE107824 — Soncin 2018
- **Current label**: Primary cytotrophoblasts (correctly labeled)
- **Article**: Isolated CTBs — "chorionic villi were minced and subjected to three sequential digestions", separated on Percoll gradient

#### GSE28551 — Sitras 2012
- **Tissue**: Central villous tissue, full-thickness minus plates
- **Article**: Biopsies from central placenta near cord insertion, fetal side, full thickness avoiding basal plate and chorionic plate
- **Verdict**: "Placenta" acceptable — full-thickness villous biopsy

#### GSE37901 — Uuskula 2012
- **Tissue**: Full-thickness placental biopsy (all layers)
- **Article**: "full-thickness blocks of 1–3 cm" including "all placental layers" fetal to maternal side
- **Verdict**: "Placenta" acceptable

#### GSE73374 — Martin 2015
- **Tissue**: Full-thickness placental biopsy
- **Article**: Full-thickness biopsies from center of placenta, avoiding periphery
- **Verdict**: "Placenta" acceptable

#### GSE55439 — Wolfe 2014
- **Tissue**: Placental disc biopsies (term); products of conception (T1)
- **Article**: Term: biopsies midway between cord insertion and edge. T1: tissue from products of conception after elective termination (likely mixed villi + surrounding tissue)
- **Verdict**: "Placenta" acceptable

#### GSE37653 — Khan 2014 (excluded)
- **Tissue**: Dissected chorionic villi
- **Article**: "human first trimester placental villi" — dissected villous portions from decidual and embryonic tissue under microscope
- **Verdict**: More accurately "Chorionic villi" than "Placenta"

#### GSE9984 — Mikheev 2008
- **Tissue**: Dissected placental villi
- **Article**: "the placental villus could be dissected"
- **Verdict**: More accurately "Chorionic villi" than "Placenta"

#### GSE22490 — Rull 2013
- **Tissue**: Trophoblastic tissue from products of conception
- **Article**: "homogenized placental tissue (containing trophoblastic and decidual material)"; figure legend calls it "Trophoblastic tissue"
- **Verdict**: "Placenta" acceptable but imprecise — contains decidual contamination

#### GSE93520 — Zhang 2017
- **Current label**: Chorionic villi (correctly labeled)
- **Article**: Chorionic villus tissue from first trimester pregnancies

### Other tissue labels verified

#### GSE70102 — Bianco 2016
- **Current label**: Basal Plate (correctly labeled)
- **Article**: "The basal plate was dissected from the placenta proper, rinsed in PBS, and diced into approximately 3x3-mm pieces"

#### GSE91189 — Garrido-Gomez 2017
- **Current label**: Smooth chorion cytotrophoblasts (correctly labeled)
- **Article**: Cytotrophoblasts isolated from smooth chorion (chorion laeve) by laser capture microdissection

### Investigated 2026-05-10

## Excluded datasets

### GSE37653 — Khan 2014
- **Excluded**: All 25 samples, status "Excluded due to poor quality"
- **Source**: Visual QC of chip images — see `one_off_scripts/pca_gse37653_good_soso.R` for quality ratings
- **Article**: `articles/references/dataset_articles/GSE37653_Khan_2014.pdf`
- **Modified**: 2026-05-10

## Changelog

### 2026-05-10
- GSE37653: excluded all 25 samples (poor chip quality)
- GSE100051: cleaned text GA to integers, cleared "Term" GA values
- GSE107824: set exact GA from GEO sample titles
- GSE122214: set bounds 7–8 from Zhao 2019 article
- GSE18044: set exact GA from GEO series matrix, set diagnosis (Healthy / Healthy, Smoking)
- GSE27272: set exact GA from GEO series matrix, set diagnosis (Healthy / Healthy, Smoking)
- GSE28551: set bounds from Sitras 2012 article (T1: 9–12, Term: 38–40)
- GSE37901: fixed corrupted GA values from article gestational days
- GSE55439: set bounds from Wolfe 2014 article (T1: 10–12, Term: 38–40)
- GSE70102: refined bounds per-group from Bianco 2016 article
- GSE7434: set diagnosis from GEO titles (Healthy / Healthy, Smoking)
- GSE91189: set Late Preterm bounds 25–34 from Garrido-Gomez 2017
- GSE9984: set bounds from Mikheev 2008 article (T1: 6–8, T2: 15–16, Term: 37–41)
- GSE18044: set 1 sample bounds 35–42 from Bruchova 2010 Table 1
- GSE12767: changed Biological.Specimen from "Chorion" to "Chorionic Villus Sampling" (CVS biopsies)
- Tissue relabel: "Chorion" → "Chorionic villi" for GDS2080, GSE10588, GSE12216, GSE13155, GSE24129, GSE30186, GSE36083, GSE44711, GSE54618, GSE57050, and 22 samples with secondaryaccession "_". All confirmed as dissected villous tissue by articles/GEO, not chorion membrane. GSE73685 kept as "Chorion" (actual chorion laeve).
- Tissue relabel: "Placenta" → "Chorionic villi" for GSE28551 (central villous tissue, plates removed per Sitras 2012), GSE9984 (dissected placental villi per Mikheev 2008), GSE37653 (dissected first trimester villi per Khan 2014, excluded)
