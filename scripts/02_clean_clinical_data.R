# 02_clean_clinical_data.R
# Clean and standardize clinical annotation fields; write to data/processed/
#
# Requires: dplyr, yaml
#
# Field mapping (GDC actual name -> standardized name used downstream):
#   bcr_patient_barcode / submitter_id -> patient_id
#   sex_at_birth                       -> sex
#   age_at_index                       -> age
#   primary_diagnosis                  -> diagnosis
#   progression_or_recurrence          -> relapse_status
#   vital_status, days_to_death, days_to_last_follow_up -> used to derive survival_status
#
# Note: sample_id (biospecimen-level) and risk_group are NOT present in GDC
# clinical data and are handled separately (see comments below).

library(dplyr)
library(yaml)

config <- yaml::read_yaml("config/config.yaml")

clean_clinical <- function(project_id, raw_dir, out_dir) {
  message("Cleaning clinical data for: ", project_id)
  
  clinical <- readRDS(file.path(raw_dir, project_id, paste0(project_id, "_clinical.rds")))
  
  clinical_clean <- clinical %>%
    transmute(
      patient_id = coalesce(bcr_patient_barcode, submitter_id),
      age = as.numeric(age_at_index),
      sex = sex_at_birth,
      diagnosis = primary_diagnosis,
      vital_status = tolower(vital_status),
      days_to_death = as.numeric(days_to_death),
      days_to_last_follow_up = as.numeric(days_to_last_follow_up),
      relapse_status = progression_or_recurrence,
      treatment_response = treatments_pharmaceutical_treatment_outcome
    ) %>%
    mutate(
      # Derived field: survival_status is not a native GDC column —
      # computed here from vital_status
      survival_status = case_when(
        vital_status == "dead" ~ "deceased",
        vital_status == "alive" ~ "living",
        TRUE ~ NA_character_
      ),
      overall_survival_days = coalesce(days_to_death, days_to_last_follow_up)
    ) %>%
    filter(!is.na(patient_id))
  
  # NOTE: sample_id and risk_group are intentionally absent here.
  # - sample_id requires a separate GDCquery against the biospecimen
  #   data category (patient-level clinical data does not include it).
  # - risk_group has no direct GDC equivalent; if required, define it
  #   downstream (e.g., cytogenetic risk categories for LAML) rather
  #   than expecting it from raw clinical data.
  
  required_fields <- c("patient_id", "vital_status", "overall_survival_days")
  missing_summary <- sapply(clinical_clean[required_fields], function(x) mean(is.na(x)))
  message("Missingness in required fields:\n",
          paste(names(missing_summary), round(missing_summary, 3), sep = ": ", collapse = "\n"))
  
  dir.create(file.path(out_dir, project_id), recursive = TRUE, showWarnings = FALSE)
  saveRDS(clinical_clean, file.path(out_dir, project_id, paste0(project_id, "_clinical_clean.rds")))
  clinical_clean
}

for (role_name in names(config$datasets)) {
  ds <- config$datasets[[role_name]]
  clean_clinical(ds$name, config$paths$data_raw, config$paths$data_processed)
}

message("Done. Cleaned clinical data saved under: ", config$paths$data_processed)