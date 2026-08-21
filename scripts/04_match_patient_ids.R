# 04_match_patient_ids.R
# Matches and merges clinical and genomic records for the PRIMARY dataset (TCGA-BRCA).
#
# Checklist covered:
#   1. Identify the shared patient/sample ID format
#   2. Standardise ID strings
#   3. Match clinical and genomic records
#   4. Report matched and unmatched records
#   5. Create a merged clinical-genomics table
#
# Expected outputs:
#   data/processed/clinical_genomics_merged.csv
#   outputs/tables/id_matching_summary.csv

library(dplyr)
library(yaml)

config <- yaml::read_yaml("config/config.yaml")

# --- 1. Load cleaned clinical and genomic tables ---
clinical <- read.csv("data/processed/clinical_cleaned.csv", stringsAsFactors = FALSE)
genomics <- read.csv("data/processed/genomics_cleaned.csv", stringsAsFactors = FALSE)

message("Loaded clinical: ", nrow(clinical), " records.")
message("Loaded genomics: ", nrow(genomics), " records.")

# --- 1b. Identify the shared ID format ---
# Both tables use TCGA patient barcodes (e.g. "TCGA-A7-A0DC") as patient_id —
# clinical is one row per patient, genomics is one row per mutation record
# (multiple rows per patient), so the join key is patient_id in both.

# --- 2. Standardise ID strings ---
# Trim whitespace and uppercase to eliminate any casing/whitespace mismatches
# before matching.
clinical <- clinical %>%
  mutate(patient_id = toupper(trimws(patient_id)))

genomics <- genomics %>%
  mutate(patient_id = toupper(trimws(patient_id)))

clinical_ids <- unique(clinical$patient_id)
genomic_ids <- unique(genomics$patient_id)

# --- 3. Match clinical and genomic records ---
matched_ids <- intersect(clinical_ids, genomic_ids)
unmatched_clinical_ids <- setdiff(clinical_ids, genomic_ids)
unmatched_genomic_ids <- setdiff(genomic_ids, clinical_ids)

# --- 4. Report matched and unmatched records ---
n_clinical <- length(clinical_ids)
n_genomic <- length(genomic_ids)
n_matched <- length(matched_ids)
n_unmatched_clinical <- length(unmatched_clinical_ids)
n_unmatched_genomic <- length(unmatched_genomic_ids)
pct_matched <- round(100 * n_matched / n_clinical, 2)

id_matching_summary <- data.frame(
  number_of_clinical_records = n_clinical,
  number_of_genomic_records = n_genomic,
  number_of_matched_records = n_matched,
  number_unmatched_in_clinical_data = n_unmatched_clinical,
  number_unmatched_in_genomic_data = n_unmatched_genomic,
  percentage_matched = pct_matched
)

message("Matched patients: ", n_matched, " (", pct_matched, "% of clinical cohort)")
message("Unmatched in clinical (no mutation data): ", n_unmatched_clinical)
message("Unmatched in genomic (no clinical record): ", n_unmatched_genomic)

# --- 5. Create a merged clinical-genomics table ---
# Inner join: keep only patients present in BOTH tables (one row per mutation,
# clinical fields repeated across each patient's mutation rows).
merged <- genomics %>%
  inner_join(clinical, by = "patient_id", suffix = c("_genomic", "_clinical"))

# --- Save outputs ---
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

write.csv(merged, "data/processed/clinical_genomics_merged.csv", row.names = FALSE)
write.csv(id_matching_summary, "outputs/tables/id_matching_summary.csv", row.names = FALSE)

message("✅ Saved: data/processed/clinical_genomics_merged.csv (", nrow(merged), " rows)")
message("✅ Saved: outputs/tables/id_matching_summary.csv")

print(id_matching_summary)
print(head(merged))