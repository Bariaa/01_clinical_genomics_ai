# 03_process_mutation_data.R
# Loads, cleans, and summarizes mutation data for the PRIMARY dataset (TCGA-BRCA).
#
# Checklist covered:
#   1. Load the mutation/genomic data
#   2. Standardise gene names
#   3. Standardise patient/sample IDs
#   4. Classify mutation types
#   5. Count mutations per patient/sample
#   6. Identify most frequently mutated genes
#   7. Save a cleaned mutation table
#
# Expected outputs:
#   data/processed/genomics_cleaned.csv
#   outputs/tables/top_mutated_genes.csv
#   outputs/tables/mutation_counts_per_patient.csv
#   outputs/tables/variant_classification_summary.csv

library(dplyr)
library(yaml)

config <- yaml::read_yaml("config/config.yaml")
primary_id <- config$datasets$primary$name   # "TCGA-BRCA"

# --- 1. Load the mutation/genomic data ---
maf_path <- file.path(config$paths$data_raw, primary_id, paste0(primary_id, "_maf.rds"))
maf_raw <- readRDS(maf_path)
maf_df <- if (inherits(maf_raw, "MAF")) maf_raw@data else as.data.frame(maf_raw)
message(primary_id, ": loaded ", nrow(maf_df), " raw mutation records.")

# --- 2. Standardise gene names ---
# Hugo_Symbol is already the standard gene naming convention; trim whitespace
# and uppercase for consistency in case of any formatting inconsistencies.
maf_df <- maf_df %>%
  mutate(gene = toupper(trimws(as.character(Hugo_Symbol))))

# --- 3. Standardise patient/sample IDs ---
# TCGA barcodes: first 12 characters = patient ID; full barcode = sample ID.
maf_df <- maf_df %>%
  mutate(
    patient_id = substr(Tumor_Sample_Barcode, 1, 12),
    sample_id = as.character(Tumor_Sample_Barcode)
  )

# --- 4. Classify mutation types ---
# Variant_Classification (e.g. Missense_Mutation, Nonsense_Mutation) and
# Variant_Type (SNP, INS, DEL) are the standard MAF classification fields.
cleaned <- maf_df %>%
  transmute(
    dataset = primary_id,
    patient_id,
    sample_id,
    gene,
    chromosome = as.character(Chromosome),
    start_position = as.character(Start_Position),
    variant_classification = as.character(Variant_Classification),
    variant_type = as.character(Variant_Type),
    reference_allele = as.character(Reference_Allele),
    tumour_allele = as.character(Tumor_Seq_Allele2),
    protein_change = if ("HGVSp_Short" %in% names(maf_df)) as.character(HGVSp_Short) else NA_character_
  )

# --- 5. Count mutations per patient/sample ---
mutation_counts_per_patient <- cleaned %>%
  count(patient_id, name = "mutation_count") %>%
  arrange(desc(mutation_count))

# --- 6. Identify most frequently mutated genes ---
top_mutated_genes <- cleaned %>%
  group_by(gene) %>%
  summarise(
    mutated_samples = n_distinct(patient_id),
    total_mutations = n()
  ) %>%
  arrange(desc(mutated_samples))

# Bonus: variant classification summary (useful QC/reporting table)
variant_classification_summary <- cleaned %>%
  count(variant_classification, name = "n_mutations") %>%
  arrange(desc(n_mutations))

# --- 7. Save all outputs ---
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

write.csv(cleaned, "data/processed/genomics_cleaned.csv", row.names = FALSE)
write.csv(top_mutated_genes, "outputs/tables/top_mutated_genes.csv", row.names = FALSE)
write.csv(mutation_counts_per_patient, "outputs/tables/mutation_counts_per_patient.csv", row.names = FALSE)
write.csv(variant_classification_summary, "outputs/tables/variant_classification_summary.csv", row.names = FALSE)

message("✅ Saved: data/processed/genomics_cleaned.csv (", nrow(cleaned), " rows)")
message("✅ Saved: outputs/tables/top_mutated_genes.csv (", nrow(top_mutated_genes), " genes)")
message("✅ Saved: outputs/tables/mutation_counts_per_patient.csv (", nrow(mutation_counts_per_patient), " patients)")
message("✅ Saved: outputs/tables/variant_classification_summary.csv")

print(head(top_mutated_genes, 10))
print(head(mutation_counts_per_patient, 10))
print(variant_classification_summary)