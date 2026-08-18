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
  ## Extracted Clinical Metadata

A sample-level clinical metadata extraction (`outputs/tables/clinical_metadata_extracted.csv`)
was generated via `scripts/extract_clinical_metadata.R`, containing the following fields:
`patient_id`, `sample_id`, `age`, `sex`, `diagnosis`, `vital_status`, `days_to_death`,
`days_to_last_follow_up`, `survival_status`, `treatment_response`, `relapse_status`, `risk_group`.

**Granularity:** one row per sample (patients with multiple samples appear multiple times).

| Dataset | Sample-level rows | % with sample_id | Notes |
|---|---|---|---|
| TCGA-BRCA | 3397 | 100% | |
| TCGA-LUAD | 1781 | 100% | |
| TCGA-LAML | 697 | 100% | |

**Field availability notes:**
- `patient_id`, `sample_id`, `age`, `sex`, `diagnosis`, `vital_status`, `days_to_death`,
  `days_to_last_follow_up`, `survival_status` — populated directly or derived from GDC clinical/biospecimen data.
- `treatment_response`, `relapse_status` — sparsely populated; most values are `N/A`
  as these are inconsistently recorded in GDC across patients.
- `risk_group` — always `N/A`; no direct GDC equivalent exists. Would need to be
  computed separately (e.g., cytogenetic risk criteria for LAML) if required downstream.
  
  ## Extracted Mutation Metadata

A row-level genomic mutation extraction (`outputs/tables/mutation_metadata_extracted.csv`)
was generated via `scripts/03b_extract_mutation_metadata.R` from the downloaded MAF
(Masked Somatic Mutation) files, containing the following fields:
`patient_id`, `sample_id`, `gene`, `chromosome`, `start_position`,
`variant_classification`, `variant_type`, `reference_allele`, `tumour_allele`,
`protein_change`, `copy_number_status`, `amplification`, `deletion`, `mutation_burden`.

**Granularity:** one row per individual mutation record (not aggregated by gene or sample).

| Dataset | Mutation rows | Samples covered | Notes |
|---|---|---|---|
| TCGA-BRCA | 89,568 | 992 | |
| TCGA-LUAD | 194,729 | 618 | Highest mutation burden — consistent with smoking-related mutagenesis typical of lung adenocarcinoma |
| TCGA-LAML | 3,900 | 153 | Lowest mutation burden — consistent with AML's typically low somatic mutation rate relative to solid tumors |

**Field availability notes:**
- `gene`, `chromosome`, `start_position`, `variant_classification`, `variant_type`,
  `reference_allele`, `tumour_allele`, `protein_change` — populated directly from MAF data.
- `mutation_burden` — computed (mutation count per sample), not a native MAF field.
- `copy_number_status`, `amplification`, `deletion` — always `N/A`. These require a
  separate Copy Number Variation query (different GDC data category from Simple
  Nucleotide Variation) and are out of scope for this extraction script. If needed,
  a follow-up script (e.g., `03c_extract_cnv_metadata.R`) would be required.