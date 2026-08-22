# 07_build_ml_feature_table.R
# Converts clinical and genomic information into a machine-learning-ready
# feature table for the PRIMARY dataset (TCGA-BRCA).
#
# Example format:
#   patient_id | age | sex | mutation_burden | TP53_mutated | KRAS_mutated | PIK3CA_mutated | outcome
#
# Expected output:
#   data/processed/ml_feature_table.csv

library(dplyr)
library(tidyr)
library(yaml)

config <- yaml::read_yaml("config/config.yaml")

clinical <- read.csv("data/processed/clinical_cleaned.csv", stringsAsFactors = FALSE)
genomics <- read.csv("data/processed/genomics_cleaned.csv", stringsAsFactors = FALSE)
mutation_counts <- read.csv("outputs/tables/mutation_counts_per_patient.csv", stringsAsFactors = FALSE)
top_genes <- read.csv("outputs/tables/top_mutated_genes.csv", stringsAsFactors = FALSE)

message("Loaded clinical: ", nrow(clinical), " | genomics: ", nrow(genomics), " rows")

# =========================================================================
# CLINICAL FEATURES
# =========================================================================
clinical_features <- clinical %>%
  transmute(
    patient_id,
    age = age_at_index,
    sex = sex_at_birth,
    outcome = vital_status,                            # alive / dead
    risk_group = NA_character_,                          # not available in GDC clinical data
    treatment_response = if ("treatments_pharmaceutical_treatment_outcome" %in% names(clinical)) {
      treatments_pharmaceutical_treatment_outcome
    } else NA_character_,
    relapse_status = if ("progression_or_recurrence" %in% names(clinical)) {
      progression_or_recurrence
    } else NA_character_
  )

# =========================================================================
# GENOMIC FEATURES
# =========================================================================

# --- Mutation burden (per patient) ---
burden <- mutation_counts %>%
  rename(mutation_burden = mutation_count)

# --- Gene mutation status: TP53, KRAS, PIK3CA (per example format) + top-10 genes ---
genes_of_interest <- unique(c("TP53", "KRAS", "PIK3CA", head(top_genes$gene, 10)))

gene_status <- genomics %>%
  filter(gene %in% genes_of_interest) %>%
  distinct(patient_id, gene) %>%
  mutate(mutated = 1, gene_col = paste0(gene, "_mutated")) %>%
  select(patient_id, gene_col, mutated) %>%
  pivot_wider(names_from = gene_col, values_from = mutated, values_fill = 0)

# --- Driver gene mutation status: known BRCA driver genes ---
driver_genes <- c("TP53", "PIK3CA", "GATA3", "CDH1", "MAP3K1", "PTEN", "AKT1", "ARID1A")

driver_status <- genomics %>%
  filter(gene %in% driver_genes) %>%
  distinct(patient_id) %>%
  mutate(driver_gene_mutated = 1)

# --- Variant class counts: wide table, one column per variant classification ---
variant_class_counts <- genomics %>%
  count(patient_id, variant_classification) %>%
  pivot_wider(
    names_from = variant_classification,
    values_from = n,
    values_fill = 0,
    names_prefix = "vc_"
  )

# --- Pathway-level mutation flags ---
pathway_genes <- list(
  pathway_pi3k_akt = c("PIK3CA", "PTEN", "AKT1"),
  pathway_tp53 = c("TP53"),
  pathway_chromatin_remodeling = c("ARID1A", "KMT2C"),
  pathway_wnt = c("CTNNB1", "APC")
)

pathway_flags <- genomics %>%
  distinct(patient_id, gene) %>%
  { pw <- .
  for (pw_name in names(pathway_genes)) {
    pw[[pw_name]] <- as.integer(pw$gene %in% pathway_genes[[pw_name]])
  }
  pw
  } %>%
  group_by(patient_id) %>%
  summarise(across(starts_with("pathway_"), max), .groups = "drop")

# =========================================================================
# COMBINE INTO FINAL ML-READY FEATURE TABLE
# =========================================================================
ml_feature_table <- clinical_features %>%
  left_join(burden, by = "patient_id") %>%
  left_join(gene_status, by = "patient_id") %>%
  left_join(driver_status, by = "patient_id") %>%
  left_join(variant_class_counts, by = "patient_id") %>%
  left_join(pathway_flags, by = "patient_id") %>%
  mutate(
    across(ends_with("_mutated"), ~ replace_na(., 0)),
    across(starts_with("vc_"), ~ replace_na(., 0)),
    across(starts_with("pathway_"), ~ replace_na(., 0)),
    driver_gene_mutated = replace_na(driver_gene_mutated, 0),
    mutation_burden = replace_na(mutation_burden, 0)
  ) %>%
  relocate(patient_id, age, sex, mutation_burden, TP53_mutated, KRAS_mutated, PIK3CA_mutated, outcome)

message("Final ML feature table: ", nrow(ml_feature_table), " patients, ",
        ncol(ml_feature_table), " features")

# =========================================================================
# SAVE OUTPUT
# =========================================================================
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

write.csv(ml_feature_table, "data/processed/ml_feature_table.csv", row.names = FALSE)

message("✅ Saved: data/processed/ml_feature_table.csv (", nrow(ml_feature_table), " x ", ncol(ml_feature_table), ")")

print(head(ml_feature_table[, c("patient_id", "age", "sex", "mutation_burden",
                                "TP53_mutated", "KRAS_mutated", "PIK3CA_mutated", "outcome")]))