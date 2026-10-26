# 13_integrate_epigenomics.R
# Integrate clinical, genomic and epigenomic layers into one patient-level dataset.
#
# Steps:
#   1. Load cleaned methylation matrix and the UUID -> patient mapping
#   2. Collapse patients with multiple aliquots, keeping the highest-quality
#      sample (lowest pre-imputation missingness)
#   3. Transpose to patient-level rows
#   4. Join clinical and genomic summary features from Month 1
#   5. Restrict to the three-way intersection
#
# Outputs:
#   data/processed/multiomics_integrated.csv   (gitignored - large)
#   outputs/tables/multiomics_integration_summary.csv

library(data.table)
library(dplyr)
library(yaml)

config <- yaml::read_yaml("config/config.yaml")
primary_id <- config$datasets$primary$name
proc_dir <- config$paths$data_processed

# --- 1. Load inputs ---
beta <- readRDS(file.path(proc_dir, paste0(primary_id, "_methylation_cleaned.rds")))
id_map <- readRDS(file.path(proc_dir, paste0(primary_id, "_methylation_id_map.rds")))
three_way <- readRDS(file.path(proc_dir, paste0(primary_id, "_three_way_patient_ids.rds")))

message("Methylation matrix: ", paste(dim(beta), collapse = " x "))
message("Three-way patients: ", length(three_way))

# --- 2. Collapse duplicate aliquots: keep lowest pre-imputation missingness ---
raw_beta <- readRDS(file.path(config$paths$data_raw, primary_id,
                              paste0(primary_id, "_methylation_betas.rds")))
sample_missing <- colMeans(is.na(raw_beta))

map <- id_map %>%
  filter(file_uuid %in% colnames(beta)) %>%
  mutate(pct_missing = 100 * sample_missing[file_uuid])

n_dup_patients <- sum(duplicated(map$patient_id))
message("Patients with multiple aliquots: ", length(unique(map$patient_id[duplicated(map$patient_id)])))

best <- map %>%
  group_by(patient_id) %>%
  slice_min(pct_missing, n = 1, with_ties = FALSE) %>%
  ungroup()

message("Aliquots dropped as lower-quality duplicates: ", nrow(map) - nrow(best))

# --- 3. Transpose to patient-level ---
beta_best <- beta[, best$file_uuid, drop = FALSE]
colnames(beta_best) <- best$patient_id

meth_patient <- as.data.frame(t(beta_best))
meth_patient$patient_id <- rownames(meth_patient)
meth_patient <- meth_patient %>% relocate(patient_id)
message("Patient-level methylation: ", paste(dim(meth_patient), collapse = " x "))

# --- 4. Clinical features (Month 1) ---
clinical <- read.csv("data/processed/clinical_cleaned.csv", stringsAsFactors = FALSE) %>%
  mutate(patient_id = toupper(trimws(patient_id)))

key_clinical <- intersect(
  c("patient_id", "age_at_index", "sex_at_birth", "primary_diagnosis",
    "vital_status", "overall_survival_days", "survival_event"),
  names(clinical)
)
clinical_small <- clinical[, key_clinical, drop = FALSE]

# --- 4b. Genomic summary features (Month 1) ---
genomics <- read.csv("data/processed/genomics_cleaned.csv", stringsAsFactors = FALSE) %>%
  mutate(patient_id = toupper(trimws(patient_id)))

genomic_summary <- genomics %>%
  group_by(patient_id) %>%
  summarise(
    mutation_burden = n(),
    TP53_mutated   = as.integer(any(gene == "TP53")),
    PIK3CA_mutated = as.integer(any(gene == "PIK3CA")),
    .groups = "drop"
  )

# --- 5. Join and restrict to three-way intersection ---
integrated <- meth_patient %>%
  inner_join(clinical_small, by = "patient_id") %>%
  inner_join(genomic_summary, by = "patient_id") %>%
  filter(patient_id %in% three_way)

message("Integrated dataset: ", paste(dim(integrated), collapse = " x "))

n_meth_features <- sum(grepl("^cg", names(integrated)))
n_clin_features <- length(key_clinical) - 1
n_gen_features  <- ncol(genomic_summary) - 1

integration_summary <- data.frame(
  n_patients_integrated = nrow(integrated),
  n_total_columns = ncol(integrated),
  n_methylation_features = n_meth_features,
  n_clinical_features = n_clin_features,
  n_genomic_features = n_gen_features,
  n_methylation_patients_before_join = nrow(meth_patient),
  n_duplicate_aliquots_dropped = nrow(map) - nrow(best),
  pct_of_clinical_cohort = round(100 * nrow(integrated) / nrow(clinical), 2)
)

# --- Save ---
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
data.table::fwrite(integrated, "data/processed/multiomics_integrated.csv")
write.csv(integration_summary,
          "outputs/tables/multiomics_integration_summary.csv", row.names = FALSE)
saveRDS(integrated, file.path(proc_dir, "multiomics_integrated.rds"))

message("\n✅ Saved: data/processed/multiomics_integrated.csv (",
        round(file.size("data/processed/multiomics_integrated.csv")/1024^2, 1), " MB)")
message("✅ Saved: outputs/tables/multiomics_integration_summary.csv")
print(t(integration_summary))
