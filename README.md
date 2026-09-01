# # Clinical Genomics AI: TCGA Multi-Cancer Pipeline Generalizability Study

Clinical genomics AI pipeline: data preparation, cohort matching, feature
engineering, and baseline modeling across TCGA datasets, with a secondary
dataset used to test **pipeline generalizability** (not model generalizability).

## Study Design

| Role | Dataset | Purpose |
|---|---|---|
| Primary (development) | TCGA-BRCA | Build and debug the pipeline against the largest, most feature-complete cohort |
| Secondary (validation) | TCGA-LUAD | Apply the unmodified pipeline to a structurally similar solid-tumor cohort |
| Tertiary (stress test) | TCGA-LAML | Apply the unmodified pipeline to a structurally distinct liquid-tumor cohort to probe robustness and failure modes |

Full justification is documented in `outputs/reports/dataset_selection_report.md`.

## Repository Structure

01_clinical_genomics_ai/
├── data/
│ ├── raw/ # Original downloads — never edited directly
│ ├── processed/ # Cleaned, analysis-ready data
│ └── example/ # Small example data for testing scripts
├── scripts/
│ ├── 01_download_or_prepare_data.R
│ ├── 02_clean_clinical_data.R
│ ├── 03_process_mutation_data.R
│ ├── 04_match_patient_ids.R
│ ├── 05_generate_summary_tables.R
│ ├── 06_generate_figures.R
│ ├── 07_build_ml_feature_table.R
│ ├── 08_train_baseline_model.R
│ └── 09_test_second_dataset.R
├── notebooks/
│ ├── 01_exploratory_data_review.Rmd
│ └── 02_model_exploration.Rmd
├── outputs/
│ ├── tables/
│ ├── figures/
│ ├── models/
│ └── reports/
│ ├── data_overview.md
│ └── dataset_selection_report.md
├── config/
│ └── config.yaml
├── environment/
│ └── environment.yml
├── renv.lock # generated after running renv::init()
└── README.md

## Input Data Requirements

This pipeline requires raw TCGA data to be downloaded before running any
processing scripts:
- **Clinical data**: `data/raw/<DATASET>/<DATASET>_clinical.rds` (via `TCGAbiolinks::GDCquery_clinic()`)
- **Mutation data (MAF)**: `data/raw/<DATASET>/<DATASET>_maf.rds` (via `TCGAbiolinks::GDCquery()` + `GDCdownload()`)

Downloads require an internet connection and, on Windows, the standalone
[GDC Data Transfer Tool](https://gdc.cancer.gov/access-data/gdc-data-transfer-tool)
for reliable large-file downloads (`method = "client"` in `GDCdownload()`).
Raw data is not committed to this repository — each user must download it
independently (see `.gitignore`).

## Installation Instructions

1. Install [R](https://cran.r-project.org/) (4.3+) and [RStudio](https://posit.co/download/rstudio-desktop/)
2. Clone this repository
3. Open `01_clinical_genomics_ai.Rproj` in RStudio
4. Install `renv` and restore the locked package environment (see Setup below)
5. On Windows, install [Rtools](https://cran.rstudio.com/bin/windows/Rtools/) (see Prerequisites below)

## Conda Setup Instructions

An alternate conda-based environment is provided in `environment/environment.yml`,
useful when R needs to run via conda (e.g. on a shared cluster) rather than a
system R install:
```bash
conda env create -f environment/environment.yml
conda activate clinical-genomics-ai
```
**Note:** `renv` (see Setup below) is the primary reproducibility method for
this project; the conda environment is a secondary option.

## Docker Setup Instructions

A first-attempt `environment/Dockerfile` is provided for containerising the
environment:
```bash
docker build -t clinical-genomics-ai -f environment/Dockerfile .
docker run -it clinical-genomics-ai
```
**Known limitation:** this Dockerfile has not yet been built or tested (Docker
was not available in the development environment used for this project) — see
the note at the top of the Dockerfile itself.

## Expected Outputs

| Script | Key outputs |
|---|---|
| 02 | `data/processed/clinical_cleaned.csv`, `outputs/tables/clinical_missingness_summary.csv` |
| 03 | `data/processed/genomics_cleaned.csv`, `outputs/tables/top_mutated_genes.csv`, `mutation_counts_per_patient.csv`, `variant_classification_summary.csv` |
| 04 | `data/processed/clinical_genomics_merged.csv`, `outputs/tables/id_matching_summary.csv` |
| 05 | `outputs/tables/mutation_burden_summary.csv`, `clinical_genomic_summary.csv` |
| 06 | `outputs/figures/*.png` (top genes, variant distribution, mutation counts, clinical grouping) |
| 07 | `data/processed/ml_feature_table.csv` |
| 08 | `outputs/tables/model_performance_summary.csv`, `feature_importance.csv`, `outputs/figures/confusion_matrix.png`, `roc_curve.png`, `outputs/models/baseline_model.rds` |
| LUAD test | `outputs/reports/luad_generalizability_test_log.csv`, `luad_generalizability_report.md` |

## Limitations

- **No strong native clinical endpoint**: `survival_status` showed near-chance
  predictive power (AUC ~0.556); `relapse_status` and `treatment_response` were
  too sparse to model. An exploratory `risk_group` endpoint (mutation-burden
  median split) was used instead — see `outputs/reports/data_overview.md`.
- **Follow-up time data gap**: `days_to_last_follow_up` is populated for only
  ~0.1% of BRCA patients via the standard clinical query — a known GDC API
  limitation, documented in `outputs/reports/data_overview.md`.
- **Dockerfile untested**: a first attempt only, not yet built/run.
- **No `--config` CLI flag**: scripts read `config/config.yaml` via a fixed path.
- **Nextflow orchestration is sequential only**, and requires WSL2 on Windows
  (Nextflow doesn't run natively on Windows) — see `workflow/nextflow/README_nextflow.md`.
- **Generalizability tested on LUAD only** (not yet extended to LAML).

## Prerequisites (Windows)
Rtools is required to build some Bioconductor dependencies from source
(e.g., `httr2`, `rlang`). Install from: https://cran.rstudio.com/bin/windows/Rtools/

## Setup

```r
install.packages("renv")
renv::init()       # detects packages used in scripts/ and creates renv.lock
```

Collaborators reproduce the exact environment with:
```r
renv::restore()
```

## Running the Pipeline

Scripts are numbered to indicate execution order:

```r
source("scripts/01_download_or_prepare_data.R")
source("scripts/02_clean_clinical_data.R")
source("scripts/03_process_mutation_data.R")
source("scripts/04_match_patient_ids.R")
source("scripts/05_generate_summary_tables.R")
source("scripts/06_generate_figures.R")
source("scripts/07_build_ml_feature_table.R")
source("scripts/08_train_baseline_model.R")
source("scripts/09_test_second_dataset.R")
```

All parameters (paths, dataset names, thresholds) are read from `config/config.yaml`
— values are not hardcoded inside individual scripts.

## Notebooks

`notebooks/` contains exploratory `.Rmd` files for interactive review — not part
of the production pipeline. Knit with:
```r
rmarkdown::render("notebooks/01_exploratory_data_review.Rmd")
```

## Data Source

TCGA data is accessed via the [GDC Data Portal](https://portal.gdc.cancer.gov/)
and the `TCGAbiolinks` R package, under open-access data use policies. Raw data
is not committed to this repository (see `.gitignore`).

## Citation

Weinstein JN, et al. *The Cancer Genome Atlas Pan-Cancer analysis project.*
Nat Genet. 2013.

## Example Run Commands

This is an R-based pipeline (not Python). Scripts currently read
`config/config.yaml` via a hardcoded path rather than a `--config`
command-line flag — a known limitation, noted below.

```bash
Rscript scripts/02_clean_clinical_data.R
Rscript scripts/03_process_mutation_data.R
Rscript scripts/04_match_patient_ids.R
Rscript scripts/05_generate_summary_tables.R
Rscript scripts/06_generate_figures.R
Rscript scripts/07_build_ml_feature_table.R
Rscript scripts/08_train_baseline_model.R
Rscript scripts/test_luad_generalizability.R   # secondary dataset test
```

**Known limitation:** scripts do not yet accept a `--config` argument;
`config/config.yaml` is read via a fixed relative path from the project
root. A future improvement would add CLI argument parsing (e.g. via the
`optparse` R package) so scripts can be pointed at alternate config files
without editing the hardcoded path.