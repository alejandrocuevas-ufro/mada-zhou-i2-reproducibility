# Data dictionary

## `published_reanalyses_2x2.csv`

| Variable | Meaning |
|---|---|
| `article` | Source publication and year |
| `group` | Diagnostic technology / subgroup used to define the meta-analysis |
| `study_arm` | Study or device/condition arm label as represented in the source extraction |
| `TP` | True positives |
| `FP` | False positives |
| `TN` | True negatives |
| `FN` | False negatives |
| `N` | Sum of TP + FP + TN + FN for that row |

For Manetas-Stavrakakis, multiple arms from the same publication are retained because this reproduces the unit structure used by the original meta-analysis. Their inclusion here reproduces the published analysis and does not constitute an independent methodological endorsement of treating those arms as independent studies.

For Pinto-Villalba, the unit reported in the source article is ultrasound examinations rather than necessarily unique patients.

## `metadta_telomerase_reconstruction.csv`

Study-level 2×2 counts used in the public `metadta` tutorial example. The accompanying script uses model parameters printed in that publication; because those parameters were rounded in print, this is primarily a check of the formula and denominator definitions rather than an exact re-fit of the original model.

## `data_provenance.csv`

Documents the source location and extraction status of each dataset used in the manuscript.
