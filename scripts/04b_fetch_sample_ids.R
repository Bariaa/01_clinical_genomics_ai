# 04b_fetch_sample_ids.R
# Fetch sample_id (biospecimen-level identifiers) via a separate GDC query.
# sample_id is NOT part of the clinical data category — it lives under
# biospecimen, and a single patient can have multiple samples (e.g., primary
# tumor + normal tissue), so this produces a sample-level table, not
# patient-level like scripts 02/07.
#
# Requires: TCGAbiolinks, dplyr, yaml

library(TCGAbiolinks)
library(dplyr)
library(yaml)

config <- yaml::read_yaml("config/config.yaml")

fetch_sample_ids <- function(project_id, out_dir) {
  message("Fetching biospecimen sample IDs for: ", project_id)
  
  biospecimen <- tryCatch(
    GDCquery_clinic(project = project_id, type = "biospecimen"),
    error = function(e) {
      message("Biospecimen query failed for ", project_id, ": ", conditionMessage(e))
      NULL
    }
  )
  
  if (is.null(biospecimen) || nrow(biospecimen) == 0) {
    message("No biospecimen data returned for ", project_id, " — skipping.")
    return(invisible(NULL))
  }
  
  sample_ids <- biospecimen %>%
    transmute(
      patient_id = coalesce(bcr_patient_barcode, submitter_id),
      sample_id = submitter_id,   # at the biospecimen level, submitter_id is the sample/aliquot barcode
      sample_type = sample_type
    ) %>%
    filter(!is.na(sample_id))
  
  message(project_id, ": ", nrow(sample_ids), " sample records found (",
          n_distinct(sample_ids$patient_id), " unique patients)")
  
  dir.create(file.path(out_dir, project_id), recursive = TRUE, showWarnings = FALSE)
  saveRDS(sample_ids, file.path(out_dir, project_id, paste0(project_id, "_sample_ids.rds")))
  sample_ids
}

for (role_name in names(config$datasets)) {
  ds <- config$datasets[[role_name]]
  fetch_sample_ids(ds$name, config$paths$data_processed)
}

message("Done. Sample ID tables saved under: ", config$paths$data_processed)