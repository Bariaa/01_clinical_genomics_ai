# Module Documentation

## What This Module Does

This module implements a clinical genomics AI pipeline for TCGA (The Cancer
Genome Atlas) data: it downloads, cleans, and merges clinical and mutation
(genomic) data for a primary cancer cohort (TCGA-BRCA), builds an ML-ready
feature table, trains and evaluates an exploratory classification model, and
tests the pipeline's generalizability on a second cancer cohort (TCGA-LUAD).
The pipeline is also wrapped as a reproducible Nextflow workflow
(`workflow/nextflow/`).

---

## What Each Script Does

### scripts/01_download_or_prepare_data.R
**Purpose:** Downloads raw clinical and mutation (MAF) data from GDC via
TCGAbiolinks for datasets defined in `config/config.yaml`.
**Outputs:** `.rds` files in `data/raw/<dataset>/`.

### scripts/02_clean_clinical_data.R
**Purpose:** Loads raw clinical data, standardises column names
(`janitor::clean_names()`), identifies the patient ID column, removes
duplicate records, formats categorical variables, prepares survival
variables, and flattens any list-type GDC fields (e.g. `sites_of_involvement`).
**Outputs:** `data/processed/clinical_cleaned.csv`,
`outputs/tables/clinical_missingness_summary.csv`.

### scripts/03_process_mutation_data.R
**Purpose:** Loads MAF mutation data, standardises gene names and patient/
sample IDs, classifies mutation types, counts mutations per patient, and
identifies the most frequently mutated genes.
**Outputs:** `data/processed/genomics_cleaned.csv`,
`outputs/tables/top_mutated_genes.csv`,
`outputs/tables/mutation_counts_per_patient.csv`,
`outputs/tables/variant_classification_summary.csv`.

### scripts/04_match_patient_ids.R
**Purpose:** Matches clinical and genomic records by patient ID, reports
matched/unmatched counts, and creates a merged clinical-genomics table
(mutation-row granularity, with a focused set of clinical columns to keep
file size manageable).
**Outputs:** `data/processed/clinical_genomics_merged.csv`,
`outputs/tables/id_matching_summary.csv`.

### scripts/05_generate_summary_tables.R
**Purpose:** Consolidates summary statistics from the cleaned/merged data:
mutation burden distribution and a combined clinical-genomic cohort overview.
**Outputs:** `outputs/tables/mutation_burden_summary.csv`,
`outputs/tables/clinical_genomic_summary.csv` (plus re-derived copies of
script 03's tables).

### scripts/06_generate_figures.R
**Purpose:** Generates visualizations: top mutated genes, variant
classification distribution, mutation count distribution per patient, and
mutation burden by clinical group (vital status).
**Outputs:** 4 PNG figures in `outputs/figures/`.

### scripts/07_build_ml_feature_table.R
**Purpose:** Converts clinical and genomic data into a single
machine-learning-ready feature table (see "How Mutation Features Were
Encoded" and "What Outcome Variable Was Used" below).
**Outputs:** `data/processed/ml_feature_table.csv`.

### scripts/08_train_baseline_model.R
**Purpose:** Trains and evaluates baseline classification models (Logistic
Regression, Random Forest; XGBoost attempted but not always stable given
sample size) against an exploratory endpoint (see "What Assumptions Were
Made" below), with a full evaluation suite.
**Outputs:** `outputs/tables/model_performance_summary.csv`,
`outputs/tables/feature_importance.csv`,
`outputs/figures/confusion_matrix.png`, `outputs/figures/roc_curve.png`,
`outputs/models/baseline_model.rds`.

### scripts/test_luad_generalizability.R
**Purpose:** Re-runs the core pipeline logic (clinical cleaning, mutation
processing, ID matching) against the secondary dataset (TCGA-LUAD) without
modification, to test whether the pipeline generalizes.
**Outputs:** `outputs/reports/luad_generalizability_test_log.csv`,
`outputs/reports/luad_generalizability_report.md`.

### workflow/nextflow/main.nf
**Purpose:** Orchestrates scripts 02-08 plus the LUAD generalizability test
as a sequential Nextflow pipeline. See
`workflow/nextflow/README_nextflow.md` for full details.

---

## What Each Output Means

| Output | Meaning |
|---|---|
| `clinical_cleaned.csv` | One row per patient, standardized clinical fields |
| `genomics_cleaned.csv` | One row per mutation record (multiple rows per patient) |
| `clinical_genomics_merged.csv` | Mutation-level rows joined with key clinical fields for matched patients |
| `top_mutated_genes.csv` | Genes ranked by number of patients carrying a mutation in them |
| `mutation_counts_per_patient.csv` | Total mutation count per patient ("mutation burden") |
| `variant_classification_summary.csv` | Counts of mutation types (missense, nonsense, silent, etc.) across the cohort |
| `id_matching_summary.csv` | How many clinical vs. genomic patient records matched/didn't match |
| `ml_feature_table.csv` | One row per patient, ready for modeling (see encoding details below) |
| `model_performance_summary.csv` | Accuracy/precision/recall/F1/ROC-AUC per model tested |
| `feature_importance.csv` | Random Forest variable importance ranking |
| `baseline_model.rds` | The best-performing trained model object (R format) |

---

## How IDs Were Matched

TCGA uses a standardized barcode format
(`TCGA-XX-XXXX-...`). The first **12 characters** of any barcode uniquely
identify the **patient**, regardless of how many additional segments follow
(sample, portion, aliquot). This project extracts `patient_id` as the first
12 characters of `Tumor_Sample_Barcode` (genomic data) or uses
`bcr_patient_barcode`/`submitter_id` directly (clinical data, already
patient-level). IDs are standardized (trimmed, uppercased) before matching
via an inner join on `patient_id`. Sample-level identifiers (`sample_id`)
are the full barcode, used where sample-level granularity matters (e.g.
biospecimen data, mutation records).

---

## How Mutation Features Were Encoded

In `scripts/07_build_ml_feature_table.R`:
- **Mutation burden**: total mutation count per patient (continuous feature).
- **Gene mutation status**: binary flags (`GENE_mutated`) for a defined set
  of genes of interest (e.g. `TP53_mutated`, `KRAS_mutated`,
  `PIK3CA_mutated`) plus the top-10 most frequently mutated genes in the
  cohort — 1 if the patient has ≥1 mutation in that gene, else 0.
- **Driver gene mutation status**: a single binary flag (`driver_gene_mutated`)
  indicating whether the patient has ≥1 mutation in a curated list of
  known BRCA driver genes (TP53, PIK3CA, GATA3, CDH1, MAP3K1, PTEN, AKT1,
  ARID1A).
- **Variant class counts**: per-patient counts of each mutation type
  (missense, nonsense, silent, etc.), one column per type.
- **Pathway-level flags**: binary indicators for curated gene-pathway
  groupings (PI3K/AKT, TP53, chromatin remodeling, Wnt).

**Important:** `patient_id` is retained only for tracking/matching and is
explicitly excluded from every model's feature set (enforced via
`setdiff()` and a `stopifnot()` check in script 08).

---

## What Outcome Variable Was Used

The dataset does not have a strong native clinical endpoint (see
Limitations). Three candidate targets were evaluated for feasibility
(`outputs/tables/prediction_target_feasibility.csv`):
- `survival_status` (vital_status: alive/dead) — feasible, but weak signal
  (AUC ~0.556).
- `relapse_status` — **infeasible**, 0/1098 usable records.
- `treatment_response` — **infeasible**, only 9/1098 usable records across
  5 categories.

**The exploratory endpoint ultimately used for the main model evaluation is
`risk_group`**: a binary label created by median-splitting `mutation_burden`
into "high_risk"/"low_risk". This is **not a validated clinical risk
classification** — it is an exploratory proxy chosen specifically because no
adequate native clinical endpoint exists in this cohort.
`mutation_burden` itself is deliberately excluded from the model's feature
set to prevent the label's own definition from leaking into the predictors
(an issue caught and fixed during development — an earlier version showed a
suspicious AUC of 1.0 before this fix).

---

## What Assumptions Were Made

- **TCGA-BRCA was selected as the primary dataset** based on largest sample
  size (1098 patients) and most complete clinical/molecular annotation
  (see