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