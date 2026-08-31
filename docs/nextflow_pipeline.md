# Nextflow Orchestration Layer

This project's R pipeline (scripts 02-08) is wrapped as a Nextflow workflow
(`main.nf`), providing a reproducible, modular execution layer on top of the
existing scripts.

## Requirements
- **WSL2** (Windows Subsystem for Linux) — Nextflow does not run natively on Windows.
- **Java** (OpenJDK 21+) — installed via `sudo apt install default-jdk`
- **Nextflow** — installed via `curl -s https://get.nextflow.io | bash`
- **R (Linux build)** with the same packages as the Windows pipeline:
  `yaml`, `dplyr`, `janitor`, `tidyr`, `ggplot2`, `caret`, `pROC`,
  `randomForest`, `TCGAbiolinks`, `maftools`.

Note: this is a **separate R installation** from the Windows/RStudio R used
during interactive development — packages must be installed independently
inside the WSL environment.

## Running the pipeline
```bash
cd /mnt/c/Cprojects/01_clinical_genomics_ai
nextflow run main.nf
```

## Design
Each process wraps one existing R script (`02_clean_clinical_data.R` through
`08_train_baseline_model.R`), executed via `Rscript` from the project root,
chained in sequence using lightweight `val` signals to enforce execution
order.

**Deliberate simplification:** scripts read/write to fixed paths (as defined
in `config/config.yaml`) rather than passing files through Nextflow's
channel system as typed artifacts. This wraps the existing, already-verified
scripts without rewriting their I/O logic.

**Future improvement path:** refactor scripts to accept explicit
input/output paths as CLI arguments, enabling proper Nextflow channel-based
file staging, caching, and parallel execution — a natural next step when
extending this pipeline to additional omics layers (e.g., copy number
variation, RNA expression), where independent branches could run in parallel
rather than strictly sequentially.

## Outputs
Running the pipeline regenerates all standard outputs
(`data/processed/`, `outputs/tables/`, `outputs/figures/`,
`outputs/models/`), plus two Nextflow-specific execution artifacts:
- `outputs/reports/nextflow_report.html` — resource usage per process
- `outputs/reports/nextflow_timeline.html` — execution timeline

## Verification
This orchestration layer was tested end-to-end: all 7 processes
(`CLEAN_CLINICAL` → `TRAIN_MODEL`) completed successfully on the primary
dataset (TCGA-BRCA), regenerating outputs identical in structure to the
original Windows/RStudio pipeline runs.