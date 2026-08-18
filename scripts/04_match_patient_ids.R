# 04_match_patient_ids.R
# Match/align patient identifiers across clinical and mutation data sources
#
# Requires: dplyr, yaml

library(dplyr)
library(yaml)

config <- yaml::read_yaml("config/config.yaml")
id_key <- config$pipeline$id_matching_key   # "patient_barcode"

match_ids <- function(project_id, processed_dir) {
  message("Matching patient IDs for: ", project_id)
  
  clinical <- readRDS(file.path(processed_dir, project_id, paste0(project_id, "_clinical_clean.rds")))
  gene_summary_path <- file.path(processed_dir, project_id, paste0(project_id, "_gene_mutation_summary.csv"))
  
  if (!file.exists(gene_summary_path)) {
    message("No mutation summary found for ", project_id, " — skipping ID matching for this dataset.")
    return(invisible(NULL))
  }
  
  matched_ids <- clinical %>% pull(!!id_key) %>% unique()
  
  message(project_id, ": ", length(matched_ids), " unique patient IDs in clinical data")
  
  saveRDS(matched_ids, file.path(processed_dir, project_id, paste0(project_id, "_matched_ids.rds")))
  matched_ids
}

for (role_name in names(config$datasets)) {
  ds <- config$datasets[[role_name]]
  match_ids(ds$name, config$paths$data_processed)
}

message("Done. Matched ID lists saved under: ", config$paths$data_processed)