# Data Overview

> Summary of the three TCGA cohorts used in this project. This should be
> regenerated/updated after running the data ingestion and cleaning scripts
> (01–02), so the numbers below reflect the actual downloaded data rather
> than placeholders.

## Cohort Summary

| Dataset | Role | N patients | Tumor type | Key annotations available |
|---|---|---|---|---|
| TCGA-BRCA | Primary | 1088 | Solid (breast) | Stage, ER/PR/HER2, histologic subtype, survival |
| TCGA-LUAD | Secondary | ~500 | Solid (lung) | Stage, smoking status, histologic subtype, survival |
| TCGA-LAML | Tertiary / stress test | ~200 | Liquid (blood/marrow) | FAB subtype, cytogenetic risk, blast %, survival |

## Field Completeness

| Field | BRCA % complete | LUAD % complete | LAML % complete |
|---|---|---|---|
| Overall survival | — | — | — |
| Vital status | — | — | — |
| Tumor stage | — | — | N/A (liquid tumor) |
| Histologic subtype | — | — | N/A (liquid tumor) |

*(Fill in the "—" values with actual percentages once `scripts/02_clean_clinical_data.R`
has been run — the completeness numbers are printed to the console when that
script runs, and can be copied here.)*

## Notes
- LAML lacks solid-tumor fields (stage, laterality, histologic subtype) by design —
  this is expected and tracked as part of the stress-test phase (see
  `dataset_selection_report.md`), not treated as missing data to impute.
- This file should be kept in sync with the actual data version used — update it
  whenever the underlying dataset is re-downloaded or re-processed.