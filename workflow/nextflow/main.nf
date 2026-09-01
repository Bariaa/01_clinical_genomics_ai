#!/usr/bin/env nextflow
/*
 * main.nf — Orchestration layer for the clinical genomics AI pipeline.
 * Location: workflow/nextflow/main.nf
 *
 * Wraps existing R scripts as sequential Nextflow processes, run against
 * the primary dataset (TCGA-BRCA), plus a generalizability test against
 * the secondary dataset (TCGA-LUAD).
 *
 * DESIGN NOTE: Each script reads/writes to fixed paths relative to the
 * project root (as defined in config/config.yaml), rather than passing
 * files through Nextflow's channel system as typed artifacts. Processes
 * are chained via simple completion signals (val emitting "done") to
 * enforce execution order. This is a deliberate simplification for this
 * first orchestration layer. See README_nextflow.md for known limitations
 * and future improvement path.
 */

params.project_dir = "${projectDir}/../.."

process CLEAN_CLINICAL_DATA {
    tag "02_clean_clinical_data"
    output:
    val "done"
    script:
    """
    cd ${params.project_dir}
    Rscript scripts/02_clean_clinical_data.R
    """
}

process PROCESS_MUTATION_DATA {
    tag "03_process_mutation_data"
    input:
    val ready
    output:
    val "done"
    script:
    """
    cd ${params.project_dir}
    Rscript scripts/03_process_mutation_data.R
    """
}

process MATCH_PATIENT_IDS {
    tag "04_match_patient_ids"
    input:
    val ready
    output:
    val "done"
    script:
    """
    cd ${params.project_dir}
    Rscript scripts/04_match_patient_ids.R
    """
}

process GENERATE_SUMMARY_TABLES {
    tag "05_generate_summary_tables"
    input:
    val ready
    output:
    val "done"
    script:
    """
    cd ${params.project_dir}
    Rscript scripts/05_generate_summary_tables.R
    """
}

process GENERATE_FIGURES {
    tag "06_generate_figures"
    input:
    val ready
    output:
    val "done"
    script:
    """
    cd ${params.project_dir}
    Rscript scripts/06_generate_figures.R
    """
}

process BUILD_ML_FEATURE_TABLE {
    tag "07_build_ml_feature_table"
    input:
    val ready
    output:
    val "done"
    script:
    """
    cd ${params.project_dir}
    Rscript scripts/07_build_ml_feature_table.R
    """
}

process TRAIN_BASELINE_MODEL {
    tag "08_train_baseline_model"
    input:
    val ready
    output:
    val "done"
    script:
    """
    cd ${params.project_dir}
    Rscript scripts/08_train_baseline_model.R
    """
}

process TEST_SECOND_DATASET {
    tag "test_luad_generalizability"
    input:
    val ready
    script:
    """
    cd ${params.project_dir}
    Rscript scripts/test_luad_generalizability.R
    """
}

workflow {
    step1 = CLEAN_CLINICAL_DATA()
    step2 = PROCESS_MUTATION_DATA(step1)
    step3 = MATCH_PATIENT_IDS(step2)
    step4 = GENERATE_SUMMARY_TABLES(step3)
    step5 = GENERATE_FIGURES(step4)
    step6 = BUILD_ML_FEATURE_TABLE(step5)
    step7 = TRAIN_BASELINE_MODEL(step6)
    TEST_SECOND_DATASET(step7)
}