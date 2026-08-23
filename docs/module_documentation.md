# Module Documentation

Detailed documentation of each script's purpose, inputs, and outputs.

## scripts/01_download_or_prepare_data.R
**Purpose:** Downloads raw clinical and mutation (MAF) data from GDC via TCGAbiolinks.
**Inputs:** Dataset names/roles from `config/config.yaml`.
**Outputs:** `.rds` files saved to `data/raw/<dataset>/`.

## scripts/02_clean_clinical_data.R
**Purpose:** Cleans and standardizes clinical fields, mapping GDC's actual
column names (e.g., `bcr_patient_barcode`, `sex_at_birth`, `age_at_index`,
`primary_diagnosis`, `progression_or_recurrence`) to standardized names used
downstream (`patient_id`, `sex`, `age`, `diagnosis`, `relapse_status`).
Also derives `survival_status` from `vital_status`, since it is not a native
GDC field.
**Inputs:** Raw clinical `.rds` files from `data/raw/`.
**Outputs:** Cleaned clinical `.rds` files in `data/processed/`.
**Note:** `sample_id` and `risk_group` are intentionally excluded — see
`scripts/04b_fetch_sample_ids.R` for sample-level IDs; `risk_group` has no
direct GDC equivalent and would need to be defined separately if required.

## scripts/03_process_mutation_data.R
**Purpose:** Processes MAF mutation data, computes gene-level mutation frequency.
**Inputs:** Raw MAF `.rds` files from `data/raw/`.
**Outputs:** Processed MAF object and gene summary CSV in `data/processed/`.

## scripts/04b_fetch_sample_ids.R
**Purpose:** Fetches sample-level identifiers (`sample_id`) via a separate
GDC biospecimen query, since sample IDs are not part of the clinical data
category and clinical data is patient-level, not sample-level.
**Inputs:** Dataset names/roles from `config/config.yaml`.
**Outputs:** Sample ID table (`patient_id`, `sample_id`, `sample_type`) saved
as `.rds` in `data/processed/<dataset>/`.
**Note:** A single patient may have multiple associated samples (e.g., tumor
and matched normal tissue), so this produces more rows than the patient-level
clinical table.

## scripts/04b_fetch_sample_ids.R
**Purpose:** Fetches sample-level identifiers (`sample_id`) via a separate
GDC biospecimen query, since sample IDs are not part of the clinical data
category and clinical data is patient-level, not sample-level.
**Inputs:** Dataset names/roles from `config/config.yaml`.
**Outputs:** Sample ID table (`patient_id`, `sample_id`, `sample_type`) saved
as `.rds` in `data/processed/<dataset>/`.
**Note:** A single patient may have multiple associated samples (e.g., tumor
and matched normal tissue), so this produces more rows than the patient-level
clinical table.

## scripts/05_generate_summary_tables.R
**Purpose:** Builds cohort summary statistics across all three datasets.
**Inputs:** Cleaned clinical data.
**Outputs:** `outputs/tables/cohort_summary.csv`.

## scripts/06_generate_figures.R
**Purpose:** Generates Kaplan-Meier survival curves per dataset.
**Inputs:** Cleaned clinical data.
**Outputs:** PNG figures in `outputs/figures/`.

## scripts/07_build_ml_feature_table.R
**Purpose:** Builds the final feature table used for baseline modeling, using
the standardized field names produced by script 02 (`age`, `sex`, `diagnosis`,
`relapse_status`, `treatment_response`) plus mutation-derived features.
**Inputs:** Cleaned clinical data (script 02 output), gene mutation summary
(script 03 output).
**Outputs:** Feature table `.rds` in `data/processed/`.

## scripts/08_train_baseline_model.R
**Purpose:** Trains a baseline logistic regression model on the primary dataset only.
**Inputs:** Feature table for TCGA-BRCA.
**Outputs:** Trained model saved to `outputs/models/baseline_model.rds`.



## scripts/08_train_baseline_model.R
**Purpose:** Trains baseline models on TCGA-BRCA against multiple candidate
prediction targets: survival_status, risk_group (proxy: mutation_burden
median split), relapse_status, treatment_response, and long_vs_short_survival.
Each target's feasibility (sample size, class balance) is checked before
training; targets below a minimum threshold (20 cases in the smaller class)
are reported but not trained. **patient_id is always excluded from the
feature set.**
**Inputs:** `data/processed/ml_feature_table.csv`, `data/processed/clinical_cleaned.csv`
**Outputs:** `outputs/tables/prediction_target_feasibility.csv`,
`outputs/tables/model_results_by_target.csv`, one `.rds` model per feasible target.

**Results:**
| Target | AUC-ROC | n |
|---|---|---|
| survival_status | 0.556 | 1097 |
| risk_group | 0.744 | 1097 |
| long_vs_short_survival | 0.635 | 151 |

**Note on risk_group:** since this target is derived from `mutation_burden`,
that predictor is deliberately excluded when training this specific model —
including it caused perfect (AUC=1.0) leakage, since the model could trivially
recover the exact rule used to define the label.

**Infeasible targets:** `relapse_status` had 0 usable records (entirely
missing/unreported in GDC clinical data for this cohort); `treatment_response`
had only 9 usable records across 5 categories — both too sparse to model.
modeling instead of binary vital status).

## scripts/09_test_second_dataset.R
**Purpose:** Applies the unmodified pipeline/model to secondary and tertiary datasets;
logs completion rate, attrition, and accuracy per dataset.
**Inputs:** Trained model, feature tables for all datasets.
**Outputs:** `outputs/reports/generalizability_test_results.csv`.