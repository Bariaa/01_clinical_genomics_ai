# 03b_extract_mutation_metadata.R
# Extracts row-level genomic mutation metadata (one row per mutation record)
# across all datasets: patient_id, sample_id, gene, chromosome, start_position,
# variant_classification, variant_type, reference_allele, tumour_allele,
# protein_change, copy_number_status, amplification, deletion, mutation_burden
#
# Note: copy_number_status, amplification, deletion are NOT available from
# MAF (Simple Nucleotide Variation) data — they require a separate Copy
# Number Variation query, which is out of scope for this script and filled
# with "N/A" accordingly. mutation_burden is computed (mutation count per
# sample), not a native MAF field.
#
# This is distinct from 03_process_mutation_data.R, which produces a
# gene-level aggregated summary rather than row-level mutation records.

library(dplyr)
library(yaml)

config <- yaml::read_yaml("config/config.yaml")

extract_mutation_data <- function(project_id, raw_dir) {
  message("Extracting mutation metadata for: ", project_id)
  
  maf_path <- file.path(raw_dir, project_id, paste0(project_id, "_maf.rds"))
  
  if (!file.exists(maf_path)) {
    message("⚠️ No MAF file found for ", project_id, " at ", maf_path, " — skipping.")
    return(NULL)
  }
  
  maf_raw <- readRDS(maf_path)
  
  # maf_raw may be a maftools MAF object or a plain data.frame depending on
  # how script 01/03 saved it — handle both
  maf_df <- if (inherits(maf_raw, "MAF")) maf_raw@data else as.data.frame(maf_raw)
  
  required_cols <- c("Tumor_Sample_Barcode", "Hugo_Symbol", "Chromosome",
                     "Start_Position", "Variant_Classification", "Variant_Type",
                     "Reference_Allele", "Tumor_Seq_Allele2")
  missing_cols <- setdiff(required_cols, names(maf_df))
  if (length(missing_cols) > 0) {
    message("⚠️ ", project_id, ": missing expected MAF columns: ",
            paste(missing_cols, collapse = ", "))
  }
  
  extracted <- maf_df %>%
    transmute(
      patient_id = substr(Tumor_Sample_Barcode, 1, 12),
      sample_id = Tumor_Sample_Barcode,
      gene = if ("Hugo_Symbol" %in% names(maf_df)) as.character(Hugo_Symbol) else NA_character_,
      chromosome = if ("Chromosome" %in% names(maf_df)) as.character(Chromosome) else NA_character_,
      start_position = if ("Start_Position" %in% names(maf_df)) as.character(Start_Position) else NA_character_,
      variant_classification = if ("Variant_Classification" %in% names(maf_df)) as.character(Variant_Classification) else NA_character_,
      variant_type = if ("Variant_Type" %in% names(maf_df)) as.character(Variant_Type) else NA_character_,
      reference_allele = if ("Reference_Allele" %in% names(maf_df)) as.character(Reference_Allele) else NA_character_,
      tumour_allele = if ("Tumor_Seq_Allele2" %in% names(maf_df)) as.character(Tumor_Seq_Allele2) else NA_character_,
      protein_change = if ("HGVSp_Short" %in% names(maf_df)) as.character(HGVSp_Short) else NA_character_,
      copy_number_status = NA_character_,
      amplification = NA_character_,
      deletion = NA_character_
    )
  
  burden <- extracted %>%
    count(sample_id, name = "mutation_burden")
  
  final <- extracted %>%
    left_join(burden, by = "sample_id") %>%
    mutate(across(everything(), as.character)) %>%
    mutate(across(everything(), ~ ifelse(is.na(.) | . %in% c("", "NA", "NaN"), "N/A", .))) %>%
    mutate(dataset = project_id, .before = patient_id)
  
  final
}

all_mutation_data <- bind_rows(
  lapply(names(config$datasets), function(role_name) {
    ds <- config$datasets[[role_name]]
    extract_mutation_data(ds$name, config$paths$data_raw)
  })
)

dir.create(config$paths$outputs_tables, recursive = TRUE, showWarnings = FALSE)
write.csv(all_mutation_data,
          file.path(config$paths$outputs_tables, "mutation_metadata_extracted.csv"),
          row.names = FALSE)

message("Done. Extracted mutation metadata saved to: ",
        file.path(config$paths$outputs_tables, "mutation_metadata_extracted.csv"))

print(head(all_mutation_data, 10))
cat("\nTotal mutation rows per dataset:\n")
print(table(all_mutation_data$dataset))