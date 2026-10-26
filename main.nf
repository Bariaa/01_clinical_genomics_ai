#!/usr/bin/env nextflow
/*
 * main.nf — Orchestration layer for the clinical genomics AI pipeline.
 *
 * Wraps existing R scripts (02-08) as sequential Nextflow processes.
 *
 * DESIGN NOTE: Each script reads/writes to fixed paths relative to the
 * project root (as defined in config/config.yaml), rather than passing
 * files through Nextflow's channel system as typed artifacts. Processes
 * are chained via simple completion signals (val emitting "done") to
 * enforce execution order, since the scripts manage their own file I/O
 * internally. This is a deliberate simplification for this first
 * orchestration layer. A future improvement would refactor scripts to
 * accept explicit input/output paths as CLI arguments, enabling proper
 * Nextflow channel-based staging and parallelization -- a natural next
 * step when extending this pipeline to additional omics layers.
 */

params.project_dir = projectDir

process CLEAN_CLINICAL {
    tag "02_clean_clinical_data"
    output:
    val "done"
    script:
    """
    cd ${params.project_dir}
    Rscript scripts/02_clean_clinical_data.R
    """
}

process PROCESS_MUTATIONS {
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

process MATCH_IDS {
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

process SUMMARY_TABLES {
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

process BUILD_FEATURES {
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

process TRAIN_MODEL {
    tag "08_train_baseline_model"
    input:
    val ready
    script:
    """
    cd ${params.project_dir}
    Rscript scripts/08_train_baseline_model.R
    """
}

workflow {
    step02 = CLEAN_CLINICAL()
    step03 = PROCESS_MUTATIONS(step02)
    step04 = MATCH_IDS(step03)
    step05 = SUMMARY_TABLES(step04)
    step06 = GENERATE_FIGURES(step05)
    step07 = BUILD_FEATURES(step06)
    TRAIN_MODEL(step07)
}