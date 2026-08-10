# 07_build_ml_feature_table.R
# Build the final feature table used for modeling; write to data/processed/
#
# Requires: dplyr, yaml

library(dplyr)
library(yaml)

config <- yaml::read_yaml("config/config.yaml")

build_feature_table <- function(project_id, processed_dir) {
  clinical <- readRDS(file.path(processed_dir, project_id, paste0(project_id, "_clinical_clean.rds")))
  
  gene_summary_path <- file.path(processed_dir, project_id, paste0(project_id, "_gene_mutation_summary.csv"))
  n_mutated_genes_top10 <- NA
  if (file.exists(gene_summary_path)) {
    gene_summary <- read.csv(gene_summary_path)
    n_mutated_genes_top10 <- nrow(head(gene_summary, 10))
  }
  
  features <- clinical %>%
    transmute(
      patient_id,
      age,
      sex,
      diagnosis,
      vital_status = as.integer(vital_status == "dead"),
      overall_survival_days,
      relapse_status,
      treatment_response,
      n_mutated_genes_top10 = n_mutated_genes_top10
    ) %>%
    filter(!is.na(overall_survival_days))
  
  saveRDS(features, file.path(processed_dir, project_id, paste0(project_id, "_feature_table.rds")))
  features
}

for (role_name in names(config$datasets)) {
  ds <- config$datasets[[role_name]]
  build_feature_table(ds$name, config$paths$data_processed)
}

message("Done. Feature tables saved under: ", config$paths$data_processed)