# Month 1 Final Report

## Project
Clinical Genomics AI Pipeline — TCGA Multi-Cancer Pipeline Generalizability Study

## Objective
Build a reproducible clinical genomics pipeline on a primary TCGA dataset
(BRCA), evaluate an exploratory classification endpoint given the absence
of a strong native clinical endpoint, and test the pipeline's
generalizability on a secondary dataset (LUAD).

## Summary of Work Completed

### 1. Dataset Selection
TCGA-BRCA selected as primary (1098 patients, most complete annotation);
TCGA-LUAD as secondary (585 patients, structurally similar solid tumor);
TCGA-LAML as tertiary stress-test dataset (200 patients, structurally
distinct liquid tumor). Full rationale in `dataset_selection_report.md`.

### 2. Data Pipeline (scripts 02-07)
- Clinical data cleaned and standardized (1098 patients, 103 columns, 0 duplicates)
- Mutation data processed (89,568 mutation records; TP53/PIK3CA top mutated
  genes — consistent with published BRCA literature)
- Clinical and genomic records matched (968/1098 patients, 88.16%)
- Summary tables and visualizations generated
- ML-ready feature table built (1098 patients × 42 features)

### 3. Exploratory Classification Model (script 08)
No strong native clinical endpoint was available: `survival_status` showed
near-chance predictive power (AUC ~0.556); `relapse_status` (0/1098 usable)
and `treatment_response` (9/1098 usable) were too sparse to model. An
exploratory `risk_group` endpoint (mutation-burden median split) was used
instead. Logistic Regression and Random Forest were trained and evaluated;
XGBoost was attempted but did not converge given sample size constraints.

**Best result:** Logistic Regression, AUC-ROC = 0.720, Accuracy = 69.4%.
A data leakage bug (mutation_burden used both as label-defining variable
and as a model feature, producing a spurious AUC of 1.0) was identified
and corrected during development.

### 4. Secondary Dataset Generalizability Test
The pipeline (scripts 02-04) was run unmodified against TCGA-LUAD. All
steps completed successfully (585 clinical records, 194,729 mutation
records, 95.21% ID match rate). Key finding: LUAD's top mutated genes
(TTN, MUC16, CSMD3...) differ substantially from BRCA's (TP53, PIK3CA),
confirming that gene-based features are dataset-specific and would need
recalculation for a LUAD-specific model. Full results in
`second_dataset_test_report.md`.

### 5. Workflow Orchestration
The full pipeline (scripts 02-08 plus the secondary dataset test) was
wrapped as a Nextflow workflow (`workflow/nextflow/main.nf`), verified with
a successful end-to-end local run (8/8 processes succeeded). Requires
WSL2 + a separate Linux R installation, documented in
`workflow/nextflow/README_nextflow.md`.

### 6. Reproducibility Infrastructure
- `config/config.yaml` — centralized parameters (paths, ID columns, outcome
  variable, feature settings, model settings)
- `environment/environment.yml` — conda environment specification
- `environment/Dockerfile` — first-attempt containerisation (untested)
- `renv.lock` — locked R package versions

## Key Findings
1. TCGA-BRCA lacks a strong native clinical endpoint suitable for
   classification; an exploratory mutation-burden-based proxy was
   necessary and is clearly documented as such, not presented as a
   validated clinical tool.
2. The pipeline's core logic (clinical cleaning, mutation processing, ID
   matching) generalizes cleanly across cancer types without code changes.
3. Gene-based and threshold-based features do **not** generalize
   automatically — they require recalculation per dataset.
4. A real data leakage issue was caught and fixed during development,
   demonstrating the pipeline evaluation process is rigorous, not just
   reporting whatever number comes out first.

## Known Limitations
- Follow-up survival time data (`days_to_last_follow_up`) is populated for
  only ~0.1% of BRCA patients via the standard GDC clinical query.
- Generalizability tested on LUAD only; LAML (tertiary stress-test dataset)
  not yet run through the pipeline.
- Dockerfile is untested (no Docker available in the development environment).
- Scripts do not yet accept a `--config` CLI argument.
- Nextflow orchestration is sequential only, using fixed file paths rather
  than typed channel-based I/O.

## Next Steps (Proposed)
- Extend generalizability testing to TCGA-LAML (tertiary stress-test dataset).
- Build dataset-specific gene/feature lists for LUAD.
- Add CLI argument parsing (`--config`) to scripts.
- Build and test the Docker container.
- Refactor Nextflow processes to use channel-based file I/O, enabling
  parallel execution and easier extension to additional omics layers.