# Clinical-Genomics Harmonisation and AI Baseline Module

## 1. Background

Publicly available cancer genomics datasets, such as those hosted by The
Cancer Genome Atlas (TCGA), combine clinical annotation with molecular
(mutation, expression, copy number) data across thousands of patients.
Harmonising these two data types — clinical records and genomic
mutation calls — into a single, analysis-ready structure is a common but
non-trivial first step in any clinical genomics AI project, complicated by
inconsistent field naming, incomplete follow-up data, and the absence of
a single obvious outcome variable to model.

This report documents Month 1 of a project to build a reproducible
clinical-genomics harmonisation pipeline, culminating in a baseline
machine learning model and an initial test of the pipeline's
generalizability to a second cancer type.

## 2. Aim

To design, implement, and validate a reproducible pipeline that:
1. Cleans and harmonises clinical and genomic data for a primary TCGA
   cancer cohort.
2. Constructs a machine-learning-ready feature table combining clinical
   and genomic information.
3. Trains and evaluates a baseline classification model, using an
   appropriately chosen outcome variable.
4. Tests whether the pipeline generalizes, unmodified, to a second cancer
   cohort.
5. Wraps the pipeline as a reproducible, orchestrated workflow (Nextflow).

## 3. Dataset Selection

Three TCGA cohorts were considered: **TCGA-BRCA** (breast cancer),
**TCGA-LUAD** (lung adenocarcinoma), and **TCGA-LAML** (acute myeloid
leukaemia). TCGA-BRCA was selected as the **primary** dataset based on its
larger sample size (1098 patients) and richer clinical/molecular
annotation. TCGA-LUAD was selected as the **secondary** dataset — a
structurally similar solid tumor from a different organ system, suitable
for testing pipeline generalizability without confounding by a
fundamentally different data schema. TCGA-LAML was designated a
**tertiary/stress-test** dataset — a hematologic malignancy with a
structurally distinct schema (no tumor staging, no solid-tumor
histology), intended to test pipeline robustness under harder,
out-of-scope conditions in later work. Full rationale is documented in
`outputs/reports/dataset_selection_report.md`.

## 4. Primary Dataset Overview

TCGA-BRCA, as retrieved via the GDC Data Portal (`TCGAbiolinks`):

| Metric | Value |
|---|---|
| Patients (clinical) | 1,098 |
| Mutation records | 89,568 |
| Patients with ≥1 mutation | 968 |
| Deceased | 13.9% (152/1098) |
| Median age | 58 |
| Top mutated genes | TP53 (341/1098 samples), PIK3CA (336/1098) |

These top-gene findings are consistent with published TCGA-BRCA mutation
landscape literature, providing biological validation that the data
extraction and processing steps are correct.

## 5. Methods

### 5.1 Clinical Data Cleaning (`scripts/02_clean_clinical_data.R`)
Raw clinical data was loaded, column names standardised
(`janitor::clean_names()`), the patient ID column identified
(`bcr_patient_barcode`/`submitter_id`), duplicate records removed (0
found), categorical variables formatted (lowercased, trimmed), survival
variables derived (`overall_survival_days`, `survival_event`), and
list-type GDC fields (e.g. `sites_of_involvement`) flattened to plain
text for CSV export. Missingness was assessed across all 103 columns.

### 5.2 Genomic Data Processing (`scripts/03_process_mutation_data.R`)
Raw MAF (Masked Somatic Mutation) files were loaded, gene names
standardised (uppercased, trimmed `Hugo_Symbol`), patient/sample IDs
extracted from `Tumor_Sample_Barcode`, mutations classified by type
(`Variant_Classification`, `Variant_Type`), mutation counts computed per
patient, and the most frequently mutated genes identified.

### 5.3 Patient/Sample ID Harmonisation (`scripts/04_match_patient_ids.R`)
TCGA barcodes follow a fixed structure in which the first 12 characters
uniquely identify the patient, regardless of downstream sample/aliquot
segments. This convention was used to standardise `patient_id` across
both clinical and genomic tables, enabling an inner join. IDs were
trimmed and uppercased prior to matching to eliminate casing/whitespace
mismatches.

### 5.4 Summary Statistics and Visualisation (`scripts/05`, `scripts/06`)
Consolidated summary tables (mutation burden distribution,
clinical-genomic cohort overview) and four visualisations (top mutated
genes, variant classification distribution, per-patient mutation count
distribution, mutation burden by vital status) were generated.

### 5.5 ML-Ready Feature Table Construction (`scripts/07_build_ml_feature_table.R`)
Clinical features (age, sex, outcome) and genomic features (mutation
burden, per-gene mutation status for genes of interest and top-10
mutated genes, a curated driver-gene flag, variant-class counts, and
curated pathway-level flags) were combined into a single patient-level
table. `patient_id` was retained solely for tracking and explicitly
excluded from all modeling steps.

### 5.6 Baseline Machine Learning Model (`scripts/08_train_baseline_model.R`)
Three candidate clinical endpoints were evaluated for feasibility
(sample size, class balance) before modeling: `survival_status`
(feasible, but weak signal), `relapse_status` (infeasible — 0/1098
usable), and `treatment_response` (infeasible — 9/1098 usable). Given
the lack of a strong native endpoint, an **exploratory classification
target, `risk_group`**, was constructed as a median split of
`mutation_burden`. `mutation_burden` itself was excluded from the model's
features to prevent label leakage — an issue identified and corrected
during development (an earlier version, which included this variable,
produced a spurious AUC of 1.0). Logistic Regression and Random Forest
were trained on an 80/20 train/test split with 5-fold cross-validation;
XGBoost was attempted but did not converge given sample size constraints.

### 5.7 Second Dataset Test (`scripts/09_test_second_dataset.R`)
The clinical cleaning, genomic processing, and ID-matching logic
(sections 5.1-5.3) were re-run, unmodified, against TCGA-LUAD to test
generalizability. Results are detailed in Section 6 and fully documented
in `outputs/reports/second_dataset_test_report.md`.

## 6. Results

| Step | BRCA (primary) | LUAD (secondary) |
|---|---|---|
| Clinical records | 1,098 | 585 |
| Clinical columns | 103 | 101 |
| Mutation records | 89,568 | 194,729 |
| ID match rate | 88.16% | 95.21% |
| Deceased rate | 13.9% | 32.1% |
| Top mutated gene | TP53 | TTN |

The pipeline's clinical-cleaning, mutation-processing, and ID-matching
logic ran **unmodified and successfully** on both datasets, confirming
the schema and barcode conventions generalize across cancer types.

## 7. Key Figures

| Figure | File | Shows |
|---|---|---|
| Top mutated genes | `top_mutated_genes_barplot.png` | Top 15 genes by mutated sample count |
| Variant classification | `variant_classification_distribution.png` | Mutation type breakdown |
| Mutations per patient | `mutations_per_patient.png` | Distribution of mutation burden |
| Mutation burden by clinical group | `mutation_burden_by_clinical_group.png` | Burden vs. vital status |
| Confusion matrix | `confusion_matrix.png` | Best model's classification results |
| ROC curve | `roc_curve.png` | Model comparison (Logistic Regression vs. Random Forest) |

## 8. Model Performance Summary

Exploratory endpoint: `risk_group` (mutation-burden median split), n=1097
(80/20 train/test split).

| Model | Accuracy | Precision | Recall | F1 | ROC-AUC |
|---|---|---|---|---|---|
| Logistic Regression | 0.694 | 0.734 | 0.622 | 0.673 | **0.720** |
| Random Forest | 0.694 | 0.750 | 0.595 | 0.663 | 0.710 |

**Top predictive features** (Random Forest importance): TP53_mutated
(dominant), driver_gene_mutated, age, PIK3CA_mutated.

## 9. Interpretation

Age, sex, and a small set of gene-mutation flags carry genuine, moderate
predictive signal (AUC ~0.72) for the exploratory mutation-burden-based
risk endpoint — meaningfully better than chance, and consistent across
two different model types (Logistic Regression, Random Forest). TP53
mutation status dominates feature importance, aligning with its
well-established role as a major prognostic marker in breast cancer. By
contrast, the same features carry almost no signal (AUC ~0.556) for
`survival_status`, the more clinically intuitive endpoint — suggesting
that in this cohort, mutation-burden-based stratification is more
learnable from available features than raw vital status, likely
influenced by BRCA's relatively low (13.9%) death rate limiting
statistical power for that target.

The generalizability test shows the pipeline's **structural** logic
(schema handling, ID matching) transfers cleanly across cancer types, but
**biological/gene-based features do not** — LUAD's mutation landscape is
substantially different (higher overall burden, different top genes),
meaning any model built on BRCA-specific gene features would need
retraining, not just re-application, for LUAD.

## 10. Limitations

- **No strong native clinical endpoint** in TCGA-BRCA; the `risk_group`
  endpoint used is an exploratory proxy, not a validated clinical risk
  classification.
- **Follow-up survival time is largely unavailable**: `days_to_last_follow_up`
  is populated for only ~0.1% of patients via the standard GDC clinical
  query, limiting true survival-time analyses to the deceased subset.
- **Gene-based features are dataset-specific** and do not transfer between
  cancer types without recalculation.
- **Generalizability tested on LUAD only**; TCGA-LAML (tertiary
  stress-test dataset) has not yet been run through the pipeline.
- **Dockerfile is a first attempt, untested** (no Docker available in the
  development environment).
- **Nextflow orchestration is sequential only**, using fixed file paths
  rather than typed channel-based I/O; scripts do not yet accept a
  `--config` command-line argument.
  
## 11. Next Steps

1. Extend the generalizability test to TCGA-LAML, the structurally
   distinct tertiary stress-test dataset.
2. Build dataset-specific gene/pathway feature lists (e.g. for LUAD:
   KRAS, EGFR, STK11, KEAP1) rather than reusing BRCA-curated lists.
3. Add CLI argument parsing (`--config`) to all scripts.
4. Build and test the Dockerfile with an actual Docker installation.
5. Refactor Nextflow processes to use channel-based file I/O, enabling
   caching, parallel execution, and easier extension to additional omics
   layers (e.g. copy number variation, RNA expression).
6. Investigate a dedicated GDC "Clinical Supplement" query to recover the
   missing follow-up time data and enable proper survival modeling.