library(TCGAbiolinks)
library(dplyr)
library(yaml)

config <- yaml::read_yaml("config/config.yaml")

extract_metadata <- function(project_id) {
  message("Extracting clinical metadata for: ", project_id)
  
  clinical <- GDCquery_clinic(project = project_id, type = "clinical")
  
  clinical_std <- clinical %>%
    transmute(
      patient_id = coalesce(bcr_patient_barcode, submitter_id),
      age = as.character(coalesce(age_at_index, age_at_diagnosis)),
      sex = if ("sex_at_birth" %in% names(clinical)) as.character(sex_at_birth) else NA_character_,
      diagnosis = if ("primary_diagnosis" %in% names(clinical)) as.character(primary_diagnosis) else NA_character_,
      vital_status = as.character(vital_status),
      days_to_death = as.character(days_to_death),
      days_to_last_follow_up = as.character(days_to_last_follow_up),
      survival_status = case_when(
        tolower(vital_status) == "dead" ~ "deceased",
        tolower(vital_status) == "alive" ~ "living",
        TRUE ~ NA_character_
      ),
      treatment_response = if ("treatments_pharmaceutical_treatment_outcome" %in% names(clinical)) {
        as.character(treatments_pharmaceutical_treatment_outcome)
      } else NA_character_,
      relapse_status = if ("progression_or_recurrence" %in% names(clinical)) {
        as.character(progression_or_recurrence)
      } else NA_character_,
      risk_group = NA_character_
    )
  
  # --- Sample-level biospecimen data (debug-friendly this time) ---
  biospecimen <- tryCatch(
    GDCquery_clinic(project = project_id, type = "biospecimen"),
    error = function(e) {
      message("❌ Biospecimen query FAILED for ", project_id, ": ", conditionMessage(e))
      NULL
    }
  )
  
  if (is.null(biospecimen)) {
    message("⚠️ ", project_id, ": biospecimen query returned NULL")
    sample_ids <- NULL
  } else if (nrow(biospecimen) == 0) {
    message("⚠️ ", project_id, ": biospecimen query returned 0 rows")
    sample_ids <- NULL
  } else {
    message("✅ ", project_id, ": biospecimen returned ", nrow(biospecimen), " rows")
    
    # sample_id is a native column here. The patient-level ID is typically
    # embedded in the first portion of submitter_id (e.g. "TCGA-A7-A0DC-01A"
    # -> patient "TCGA-A7-A0DC"), since a direct bcr_patient_barcode column
    # isn't present in this biospecimen table for these projects.
    sample_ids <- biospecimen %>%
      transmute(
        patient_id = substr(submitter_id, 1, 12),   # TCGA barcodes: first 12 chars = patient ID
        sample_id = sample_id
      ) %>%
      filter(!is.na(sample_id))
  }
  
  if (!is.null(sample_ids) && nrow(sample_ids) > 0) {
    merged <- clinical_std %>%
      left_join(sample_ids, by = "patient_id", relationship = "many-to-many")
  } else {
    merged <- clinical_std %>% mutate(sample_id = NA_character_)
  }
  
  final <- merged %>%
    select(
      patient_id, sample_id, age, sex, diagnosis, vital_status,
      days_to_death, days_to_last_follow_up, survival_status,
      treatment_response, relapse_status, risk_group
    ) %>%
    # Fix: catch true NA, empty string, AND the literal text "NA"/"NaN"
    mutate(across(everything(), as.character)) %>%
    mutate(across(everything(), ~ ifelse(is.na(.) | . %in% c("", "NA", "NaN"), "N/A", .))) %>%
    mutate(dataset = project_id, .before = patient_id)
  
  final
}

all_metadata <- bind_rows(
  extract_metadata("TCGA-BRCA"),
  extract_metadata("TCGA-LUAD"),
  extract_metadata("TCGA-LAML")
)

dir.create(config$paths$outputs_tables, recursive = TRUE, showWarnings = FALSE)
write.csv(all_metadata,
          file.path(config$paths$outputs_tables, "clinical_metadata_extracted.csv"),
          row.names = FALSE)

message("Done. Extracted metadata saved to: ",
        file.path(config$paths$outputs_tables, "clinical_metadata_extracted.csv"))

print(head(all_metadata, 10))
cat("\nTotal rows per dataset:\n")
print(table(all_metadata$dataset))