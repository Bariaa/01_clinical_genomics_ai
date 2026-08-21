# 05_generate_summary_tables.R
# Generates consolidated summary tables from the cleaned/merged data produced
# by scripts 02-04, for the PRIMARY dataset (TCGA-BRCA).
#
# Expected outputs:
#   outputs/tables/top_mutated_genes.csv
#   outputs/tables/variant_classification_summary.csv
#   outputs/tables/mutation_counts_per_patient.csv
#   outputs/tables/mutation_burden_summary.csv
#   outputs/tables/clinical_genomic_summary.csv

library(dplyr)
library(yaml)

config <- yaml::read_yaml("config/config.yaml")

genomics <- read.csv("data/processed/genomics_cleaned.csv", stringsAsFactors = FALSE)
clinical <- read.csv("data/processed/clinical_cleaned.csv", stringsAsFactors = FALSE)
merged <- read.csv("data/processed/clinical_genomics_merged.csv", stringsAsFactors = FALSE)

message("Loaded genomics: ", nrow(genomics), " rows | clinical: ", nrow(clinical),
        " rows | merged: ", nrow(merged), " rows")

# --- 1. Top mutated genes ---
top_mutated_genes <- genomics %>%
  group_by(gene) %>%
  summarise(
    mutated_samples = n_distinct(patient_id),
    total_mutations = n()
  ) %>%
  arrange(desc(mutated_samples))

# --- 2. Variant classification summary ---
variant_classification_summary <- genomics %>%
  count(variant_classification, name = "n_mutations") %>%
  arrange(desc(n_mutations))

# --- 3. Mutation counts per patient ---
mutation_counts_per_patient <- genomics %>%
  count(patient_id, name = "mutation_count") %>%
  arrange(desc(mutation_count))

# --- 4. Mutation burden summary (distribution stats across the cohort) ---
mutation_burden_summary <- mutation_counts_per_patient %>%
  summarise(
    n_patients_with_mutations = n(),
    min_mutations = min(mutation_count),
    q1_mutations = quantile(mutation_count, 0.25),
    median_mutations = median(mutation_count),
    mean_mutations = round(mean(mutation_count), 2),
    q3_mutations = quantile(mutation_count, 0.75),
    max_mutations = max(mutation_count),
    sd_mutations = round(sd(mutation_count), 2)
  )

# --- 5. Clinical-genomic summary (cohort-level overview combining both) ---
clinical_genomic_summary <- data.frame(
  dataset = config$datasets$primary$name,
  n_clinical_patients = nrow(clinical),
  n_patients_with_genomic_data = n_distinct(genomics$patient_id),
  pct_deceased = round(100 * mean(clinical$vital_status == "dead", na.rm = TRUE), 2),
  median_age = median(clinical$age_at_index, na.rm = TRUE),
  median_mutation_burden = median(mutation_counts_per_patient$mutation_count),
  top_mutated_gene = top_mutated_genes$gene[1],
  top_gene_mutated_samples = top_mutated_genes$mutated_samples[1]
)

# --- Save all outputs ---
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

write.csv(top_mutated_genes, "outputs/tables/top_mutated_genes.csv", row.names = FALSE)
write.csv(variant_classification_summary, "outputs/tables/variant_classification_summary.csv", row.names = FALSE)
write.csv(mutation_counts_per_patient, "outputs/tables/mutation_counts_per_patient.csv", row.names = FALSE)
write.csv(mutation_burden_summary, "outputs/tables/mutation_burden_summary.csv", row.names = FALSE)
write.csv(clinical_genomic_summary, "outputs/tables/clinical_genomic_summary.csv", row.names = FALSE)

message("✅ Saved all 5 summary tables to outputs/tables/")
print(mutation_burden_summary)
print(clinical_genomic_summary)