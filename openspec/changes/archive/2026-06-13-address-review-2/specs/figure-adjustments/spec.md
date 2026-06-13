## ADDED Requirements

### Requirement: Figure 2 subplot titles enlarged
The R script generating Figure 2 (`one_off_scripts/pca_before_after_combat.R`) SHALL increase subplot title font size so that the "Before ComBat — by Dataset" subtitles are clearly readable at rendered article size. The `plot.title` element_text size SHALL be increased from 15 to at least 18.

#### Scenario: Subplot titles readable at article scale
- **WHEN** Figure 2 is rendered in the article PDF
- **THEN** the subplot titles ("Before ComBat — by Dataset", "After ComBat — by Dataset", etc.) are clearly readable without zooming

### Requirement: Figure 4 correlation rounded to 3 decimal places
The R script generating Figure 4 (`scripts/integrative_analysis/article_validation/combat_sensitivity.R`) SHALL round the correlation coefficient in the scatter plot title from r=0.999134 to r=0.999, matching the precision used in the article text.

#### Scenario: Scatter plot title shows rounded correlation
- **WHEN** Figure 4 is generated
- **THEN** the scatter plot title displays "r=0.999" (not "r=0.999134")
