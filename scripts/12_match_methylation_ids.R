# 12_match_methylation_ids.R
# Harmonise methylation sample IDs with clinical and genomic patient IDs.
#
# The beta-value files are named by file UUID, not TCGA barcode, so the
# GDC query results are used to map UUID -> aliquot barcode -> patient ID
# (first 12 characters, per the TCGA barcode convention used in Month 1).
#
# Output: outputs/tables/methylation_id_matching_summary.csv

library(TCGAbiolinks)
library(dplyr)
library(yaml)

config <- yaml::read_yaml("config/config.yaml")
primary_id <- config$datasets$primary$name

# --- 1. Load cleaned methylation matrix (columns = file UUIDs) ---
beta <- readRDS(file.path(config$paths$data_processed,
                          paste0(primary_id, "_methylation_cleaned.rds")))
meth_uuids <- colnames(beta)
message("Methylation samples: ", length(meth_uuids))

# --- 2. Query GDC for the UUID -> barcode mapping ---
message("Retrieving UUID -> barcode mapping from GDC...")
query_meth <- GDCquery(
  project = primary_id,
  data.category = "DNA Methylation",
  data.type = "Methylation Beta Value",
  platform = "Illumina Human Methylation 450"
)
res <- getResults(query_meth)

# GDC returns one row per file; 'id' is the file UUID, 'cases' the aliquot barcode
mapping <- res %>%
  transmute(
    file_uuid = id,
    aliquot_barcode = cases,
    patient_id = toupper(trimws(substr(cases, 1, 12))),
    sample_type = if ("sample_type" %in% names(res)) sample_type else NA_character_
  )

matched_map <- mapping %>% filter(file_uuid %in% meth_uuids)
message("UUIDs resolved to barcodes: ", nrow(matched_map), " / ", length(meth_uuids))

# --- 3. Restrict to primary tumour samples where sample_type is available ---
# TCGA includes matched normals; mixing them with tumours would confound
# any downstream analysis, so normals are excluded here.
if (!all(is.na(matched_map$sample_type))) {
  message("\nSample type distribution:")
  print(table(matched_map$sample_type, useNA = "ifany"))
  tumour_map <- matched_map %>% filter(grepl("Primary Tumor", sample_type, ignore.case = TRUE))
  message("Primary tumour samples retained: ", nrow(tumour_map))
} else {
  message("sample_type unavailable — retaining all samples.")
  tumour_map <- matched_map
}

# --- 4. Load Month 1 clinical and genomic IDs ---
clinical <- read.csv("data/processed/clinical_cleaned.csv", stringsAsFactors = FALSE)
genomics <- read.csv("data/processed/genomics_cleaned.csv", stringsAsFactors = FALSE)

clinical_ids <- unique(toupper(trimws(clinical$patient_id)))
genomic_ids  <- unique(toupper(trimws(genomics$patient_id)))
meth_ids     <- unique(tumour_map$patient_id)

# --- 5. Overlaps ---
clin_meth      <- intersect(clinical_ids, meth_ids)
gen_meth       <- intersect(genomic_ids, meth_ids)
three_way      <- Reduce(intersect, list(clinical_ids, genomic_ids, meth_ids))

message("\n--- ID overlap ---")
message("Clinical patients:            ", length(clinical_ids))
message("Genomic patients:             ", length(genomic_ids))
message("Methylation patients:         ", length(meth_ids))
message("Clinical ∩ Methylation:       ", length(clin_meth))
message("Genomic ∩ Methylation:        ", length(gen_meth))
message("Three-way intersection:       ", length(three_way))

# --- 6. Summary table ---
id_matching_summary <- data.frame(
  n_methylation_files            = length(meth_uuids),
  n_uuids_resolved               = nrow(matched_map),
  n_primary_tumour_samples       = nrow(tumour_map),
  n_methylation_patients         = length(meth_ids),
  n_clinical_patients            = length(clinical_ids),
  n_genomic_patients             = length(genomic_ids),
  n_clinical_methylation_matched = length(clin_meth),
  n_genomic_methylation_matched  = length(gen_meth),
  n_three_way_matched            = length(three_way),
  pct_methylation_matched_to_clinical =
    round(100 * length(clin_meth) / length(meth_ids), 2),
  pct_three_way_of_clinical =
    round(100 * length(three_way) / length(clinical_ids), 2)
)

dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(id_matching_summary,
          "outputs/tables/methylation_id_matching_summary.csv", row.names = FALSE)

# Save the mapping for script 13 (gitignored alongside other processed data)
saveRDS(tumour_map, file.path(config$paths$data_processed,
                              paste0(primary_id, "_methylation_id_map.rds")))
saveRDS(three_way, file.path(config$paths$data_processed,
                             paste0(primary_id, "_three_way_patient_ids.rds")))

message("\n✅ Saved: outputs/tables/methylation_id_matching_summary.csv")
print(t(id_matching_summary))