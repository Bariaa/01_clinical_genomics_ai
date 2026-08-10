# Data Overview

> Summary of the three TCGA cohorts used in this project. Sample sizes and
> data availability below were verified directly via `TCGAbiolinks::GDCquery_clinic()`
> and `GDCquery()` against the live GDC database.

## Cohort Summary

| Dataset | Role | N patients | Tumor type | Clinical data | Mutation data | Survival data |
|---|---|---|---|---|---|---|
| TCGA-BRCA | Primary | 1098 | Solid (breast) | Yes | Yes | Yes |
| TCGA-LUAD | Secondary | 585 | Solid (lung) | Yes | Yes | Yes |
| TCGA-LAML | Tertiary / stress test | 200 | Liquid (blood/marrow) | Yes | Yes | Yes |

## Dataset Comparison Table

| Cohort | Cancer type | Clinical data | Mutation data | Survival data | Sample size | Suitable as primary dataset? | Notes |
|---|---|---|---|---|---|---|---|
| TCGA-BRCA | Breast cancer | Yes | Yes | Yes | 1098 | Yes | Largest cohort with the most complete clinical, histologic, and molecular annotation (ER/PR/HER2 status, staging, subtype). Best choice for primary/development dataset. |
| TCGA-LUAD | Lung adenocarcinoma | Yes | Yes | Yes | 585 | No — suitable as secondary/validation | Structurally similar solid tumor (staging, histology framework comparable to BRCA), independent cohort from a different organ system. Good test of generalizability across organ systems. |
| TCGA-LAML | Acute myeloid leukaemia | Yes | Yes | Yes | 200 | No — suitable as tertiary/stress test only | Hematologic malignancy — liquid tumor, structurally distinct from BRCA/LUAD (no staging, no laterality, no solid-tumor histology). Useful for testing pipeline robustness/failure modes, not standard validation. |

## Field Completeness

| Field | BRCA % complete | LUAD % complete | LAML % complete |
|---|---|---|---|
| Overall survival | — | — | — |
| Vital status | — | — | — |
| Tumor stage | — | — | N/A (liquid tumor) |
| Histologic subtype | — | — | N/A (liquid tumor) |

*(Fill in the "—" values with actual percentages once `scripts/02_clean_clinical_data.R`
has been run on the downloaded data.)*

## Notes
- LAML lacks solid-tumor fields (stage, laterality, histologic subtype) by design —
  this is expected and tracked as part of the stress-test phase (see
  `dataset_selection_report.md`), not treated as missing data to impute.
- Sample sizes verified via live GDC query on the date this file was last updated.
  Re-run the verification script periodically, as GDC cohort sizes can change
  slightly over time (e.g., embargoed cases being released).