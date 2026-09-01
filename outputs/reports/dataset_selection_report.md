# Nextflow Pipeline Run Report

## Run summary
- **Date:** 2026-08-31
- **Dataset:** TCGA-BRCA (primary)
- **Command:** `nextflow run main.nf` (from `workflow/nextflow/`)
- **Result:** ✅ All 8 processes completed successfully
- **Duration:** 3m 14s
- **CPU hours:** 0.1

## Process results
| Process | Status |
|---|---|
| CLEAN_CLINICAL_DATA | ✅ Success |
| PROCESS_MUTATION_DATA | ✅ Success |
| MATCH_PATIENT_IDS | ✅ Success |
| GENERATE_SUMMARY_TABLES | ✅ Success |
| GENERATE_FIGURES | ✅ Success |
| BUILD_ML_FEATURE_TABLE | ✅ Success |
| TRAIN_BASELINE_MODEL | ✅ Success |
| TEST_SECOND_DATASET | ✅ Success (see `outputs/reports/luad_generalizability_report.md`) |

## Environment
Run via WSL2 (Ubuntu 24.04), separate Linux R installation (R 4.3.3),
Java OpenJDK 21, Nextflow 26.04.6.

## Known limitations
1. Raw data download (GDC/`gdc-client`) is not wrapped as a Nextflow
   process — performed manually beforehand due to long runtimes and
   intermittent server timeouts requiring manual retry.
2. `TEST_SECOND_DATASET` independently re-derives LUAD's cleaned tables
   rather than reusing prior process outputs via Nextflow channels —
   scripts still use fixed project-relative I/O rather than typed channel
   artifacts (see README_nextflow.md design notes).
3. Execution is fully sequential; no parallelization implemented yet.

## Detailed execution metrics
See `outputs/reports/nextflow_report.html` and
`outputs/reports/nextflow_timeline.html` for full per-process resource
usage and timing.