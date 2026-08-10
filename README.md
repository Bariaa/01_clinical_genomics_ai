# 01_clinical_genomics_ai

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