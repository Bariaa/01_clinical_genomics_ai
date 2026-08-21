# 02_clean_clinical_data.R
# Loads raw clinical data for the PRIMARY dataset (TCGA-BRCA) and produces:
#   - data/processed/clinical_cleaned.csv
#   - outputs/tables/clinical_missingness_summary.csv
#
# Checklist covered:
#   1. Load the clinical file
#   2. Standardise column names
#   3. Identify patient/sample ID columns
#   4. Remove duplicate records
#   5. Check missing values
#   6. Format categorical variables
#   7. Prepare survival variables
#   8. Save cleaned clinical table

library(dplyr)
library(yaml)
library(janitor)

config <- yaml::read_yaml("config/config.yaml")
primary_id <- config$datasets$primary$name   # "TCGA-BRCA"

path <- file.path(config$paths$data_raw, primary_id, paste0(primary_id, "_clinical.rds"))

# --- 1. Load the clinical file ---
raw <- readRDS(path)
message(primary_id, ": loaded ", nrow(raw), " raw records.")

# --- 2. Standardise column names ---
raw <- janitor::clean_names(raw)

# --- 3. Identify patient/sample ID columns ---
id_col <- intersect(c("bcr_patient_barcode", "submitter_id"), names(raw))[1]
if (is.na(id_col)) stop("No recognizable patient ID column found for ", primary_id)

cleaned <- raw %>%
  mutate(patient_id = .data[[id_col]], dataset = primary_id) %>%
  relocate(dataset, patient_id)

# --- 4. Remove duplicate records ---
n_before <- nrow(cleaned)
cleaned <- cleaned %>% distinct(patient_id, .keep_all = TRUE)
n_removed <- n_before - nrow(cleaned)
message(primary_id, ": removed ", n_removed, " duplicate patient_id rows.")

# --- 6. Format categorical variables ---
categorical_cols <- intersect(
  c("vital_status", "sex_at_birth", "race", "ethnicity", "primary_diagnosis"),
  names(cleaned)
)
cleaned <- cleaned %>%
  mutate(across(all_of(categorical_cols), ~ trimws(tolower(as.character(.)))))

# --- Flatten any list-type columns (e.g. sites_of_involvement) before saving as CSV ---
cleaned <- cleaned %>%
  mutate(across(where(is.list), ~ sapply(., function(x) paste(x, collapse = "; "))))

# --- 7. Prepare survival variables ---
cleaned <- cleaned %>%
  mutate(
    overall_survival_days = coalesce(
      as.numeric(days_to_death),
      as.numeric(days_to_last_follow_up)
    ),
    survival_event = case_when(
      vital_status == "dead" ~ 1L,
      vital_status == "alive" ~ 0L,
      TRUE ~ NA_integer_
    )
  )

# --- 5. Check missing values ---
missingness_summary <- cleaned %>%
  summarise(across(everything(), ~ mean(is.na(.)))) %>%
  tidyr::pivot_longer(everything(), names_to = "column", values_to = "pct_missing") %>%
  arrange(desc(pct_missing))

# --- 8. Save outputs ---
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

write.csv(cleaned, "data/processed/clinical_cleaned.csv", row.names = FALSE)
write.csv(missingness_summary, "outputs/tables/clinical_missingness_summary.csv", row.names = FALSE)

message("✅ Saved: data/processed/clinical_cleaned.csv (", nrow(cleaned), " rows)")
message("✅ Saved: outputs/tables/clinical_missingness_summary.csv")

print(head(cleaned))
print(head(missingness_summary, 10))