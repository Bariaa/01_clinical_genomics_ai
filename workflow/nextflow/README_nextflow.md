# Nextflow Orchestration Layer

Wraps the clinical genomics AI pipeline (scripts in `/scripts`) as a
reproducible Nextflow workflow, including a generalizability test against
the secondary dataset (TCGA-LUAD).

## Process order
1. `CLEAN_CLINICAL_DATA`
2. `PROCESS_MUTATION_DATA`
3. `MATCH_PATIENT_IDS`
4. `GENERATE_SUMMARY_TABLES`
5. `GENERATE_FIGURES`
6. `BUILD_ML_FEATURE_TABLE`
7. `TRAIN_BASELINE_MODEL`
8. `TEST_SECOND_DATASET`

## Requirements
- WSL2 (Nextflow does not run natively on Windows)
- Java (OpenJDK 21+)
- Nextflow (`curl -s https://get.nextflow.io | bash`)
- R (Linux build, separate from Windows/RStudio R) with: yaml, dplyr,
  janitor, tidyr, ggplot2, caret, pROC, randomForest, TCGAbiolinks, maftools

## Running the pipeline
```bash
cd /mnt/c/Cprojects/01_clinical_genomics_ai/workflow/nextflow
nextflow run main.nf
```

## Design and limitations
- Scripts use fixed project-relative paths (per `config/config.yaml`)
  rather than Nextflow channel-based file I/O — a deliberate simplification
  wrapping already-verified scripts without rewriting their internals.
- `TEST_SECOND_DATASET` currently re-derives LUAD's cleaned clinical/genomic
  tables independently within its own test script, rather than reusing
  BRCA's exact process outputs through the channel system — this is a
  known limitation, not a fully generalized multi-dataset workflow yet.
- No steps are currently parallelized; execution is strictly sequential.
- **Not yet automated:** raw data download (GDC downloads via `gdc-client`)
  is a manual, separate step performed before running this workflow — it
  is not wrapped as a Nextflow process due to long runtimes, intermittent
  GDC server timeouts requiring manual retry, and the standalone
  `gdc-client` executable dependency (documented in project history).

## Future improvement path
Refactor scripts to accept explicit input/output paths as CLI arguments,
enabling proper Nextflow channel-based staging, caching, and parallel
execution — useful when extending to additional omics layers (e.g. CNV,
RNA expression) as independent parallel branches.