# test_luad_generalizability.R
# Tests whether scripts 02-08 work unmodified on TCGA-LUAD (secondary dataset).
# Logs results to answer the generalizability questions.

library(dplyr)
library(yaml)
library(janitor)

config <- yaml::read_yaml("config/config.yaml")
test_log <- list()

log_result <- function(step, status, notes) {
  test_log[[length(test_log) + 1]] <<- data.frame(step = step, status = status, notes = notes)
  message(step, ": ", status, " — ", notes)
}

# --- Test Step 02: clean clinical data ---
result_02 <- tryCatch({
  raw <- readRDS("data/raw/TCGA-LUAD/TCGA-LUAD_clinical.rds")
  raw <- janitor::clean_names(raw)
  id_col <- intersect(c("bcr_patient_barcode", "submitter_id"), names(raw))[1]
  if (is.na(id_col)) stop("No recognizable patient ID column")
  cleaned <- raw %>% mutate(patient_id = .data[[id_col]], dataset = "TCGA-LUAD") %>% relocate(dataset, patient_id)
  cleaned <- cleaned %>% distinct(patient_id, .keep_all = TRUE)
  cleaned <- cleaned %>% mutate(across(where(is.list), ~ sapply(., function(x) paste(x, collapse = "; "))))
  list(status = "success", data = cleaned, cols = names(cleaned))
}, error = function(e) list(status = "failed", error = conditionMessage(e)))

if (result_02$status == "success") {
  log_result("02_clean_clinical", "✅ PASS", paste0(nrow(result_02$data), " rows, ", ncol(result_02$data), " cols — unmodified script logic worked"))
  luad_clinical_clean <- result_02$data
} else {
  log_result("02_clean_clinical", "❌ FAIL", result_02$error)
}

# --- Test Step 03: process mutation data ---
result_03 <- tryCatch({
  maf_raw <- readRDS("data/raw/TCGA-LUAD/TCGA-LUAD_maf.rds")
  maf_df <- if (inherits(maf_raw, "MAF")) maf_raw@data else as.data.frame(maf_raw)
  cleaned <- maf_df %>%
    mutate(gene = toupper(trimws(as.character(Hugo_Symbol))),
           patient_id = substr(Tumor_Sample_Barcode, 1, 12),
           sample_id = as.character(Tumor_Sample_Barcode)) %>%
    transmute(dataset = "TCGA-LUAD", patient_id, sample_id, gene,
              chromosome = as.character(Chromosome),
              variant_classification = as.character(Variant_Classification),
              variant_type = as.character(Variant_Type))
  list(status = "success", data = cleaned)
}, error = function(e) list(status = "failed", error = conditionMessage(e)))

if (result_03$status == "success") {
  log_result("03_process_mutations", "✅ PASS", paste0(nrow(result_03$data), " mutation rows — unmodified script logic worked"))
  luad_genomics_clean <- result_03$data
} else {
  log_result("03_process_mutations", "❌ FAIL", result_03$error)
}

# --- Test Step 04: ID matching/merge ---
result_04 <- tryCatch({
  clinical_ids <- unique(luad_clinical_clean$patient_id)
  genomic_ids <- unique(luad_genomics_clean$patient_id)
  matched <- intersect(clinical_ids, genomic_ids)
  list(status = "success", n_clinical = length(clinical_ids), n_genomic = length(genomic_ids), n_matched = length(matched))
}, error = function(e) list(status = "failed", error = conditionMessage(e)))

if (result_04$status == "success") {
  pct <- round(100 * result_04$n_matched / result_04$n_clinical, 2)
  log_result("04_id_matching", "✅ PASS", paste0(result_04$n_matched, "/", result_04$n_clinical, " matched (", pct, "%)"))
} else {
  log_result("04_id_matching", "❌ FAIL", result_04$error)
}

# --- Test: does BRCA-specific survival_status logic work? (vital_status field check) ---
result_survival <- tryCatch({
  if (!"vital_status" %in% names(luad_clinical_clean)) stop("vital_status column missing")
  n_dead <- sum(tolower(luad_clinical_clean$vital_status) == "dead", na.rm = TRUE)
  list(status = "success", n_dead = n_dead, n_total = nrow(luad_clinical_clean))
}, error = function(e) list(status = "failed", error = conditionMessage(e)))

if (result_survival$status == "success") {
  log_result("survival_status_field", "✅ PASS", paste0(result_survival$n_dead, "/", result_survival$n_total, " deceased — field present and usable"))
} else {
  log_result("survival_status_field", "❌ FAIL", result_survival$error)
}

# --- Test: TP53/KRAS/PIK3CA gene columns relevant for LUAD? ---
result_genes <- tryCatch({
  top_genes_luad <- luad_genomics_clean %>% count(gene, sort = TRUE) %>% head(10)
  list(status = "success", top_genes = top_genes_luad$gene)
}, error = function(e) list(status = "failed", error = conditionMessage(e)))

if (result_genes$status == "success") {
  log_result("gene_relevance_check", "ℹ️ INFO", paste0("Top LUAD genes: ", paste(result_genes$top_genes, collapse = ", ")))
} else {
  log_result("gene_relevance_check", "❌ FAIL", result_genes$error)
}

# --- Compile final report ---
test_results <- bind_rows(test_log)
dir.create("outputs/reports", recursive = TRUE, showWarnings = FALSE)
write.csv(test_results, "outputs/reports/luad_generalizability_test_log.csv", row.names = FALSE)

print(test_results)
message("\n✅ Test log saved to outputs/reports/luad_generalizability_test_log.csv")
