# 06_generate_figures.R
# Generates visualizations from the cleaned/summarized data for the
# PRIMARY dataset (TCGA-BRCA).
#
# Expected outputs:
#   outputs/figures/top_mutated_genes_barplot.png
#   outputs/figures/variant_classification_distribution.png
#   outputs/figures/mutations_per_patient.png
#   outputs/figures/mutation_burden_by_clinical_group.png

library(dplyr)
library(ggplot2)
library(yaml)

config <- yaml::read_yaml("config/config.yaml")
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)

# --- Load required tables ---
top_mutated_genes <- read.csv("outputs/tables/top_mutated_genes.csv")
variant_classification_summary <- read.csv("outputs/tables/variant_classification_summary.csv")
mutation_counts_per_patient <- read.csv("outputs/tables/mutation_counts_per_patient.csv")
merged <- read.csv("data/processed/clinical_genomics_merged.csv", stringsAsFactors = FALSE)

# --- 1. Top mutated genes barplot ---
p1 <- top_mutated_genes %>%
  arrange(desc(mutated_samples)) %>%
  slice_head(n = 15) %>%
  ggplot(aes(x = reorder(gene, mutated_samples), y = mutated_samples)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 15 Mutated Genes — TCGA-BRCA",
       x = "Gene", y = "Number of Mutated Samples") +
  theme_minimal()

ggsave("outputs/figures/top_mutated_genes_barplot.png", p1, width = 8, height = 6)

# --- 2. Variant classification distribution ---
p2 <- variant_classification_summary %>%
  arrange(desc(n_mutations)) %>%
  ggplot(aes(x = reorder(variant_classification, n_mutations), y = n_mutations)) +
  geom_col(fill = "darkorange") +
  coord_flip() +
  labs(title = "Variant Classification Distribution — TCGA-BRCA",
       x = "Variant Classification", y = "Number of Mutations") +
  theme_minimal()

ggsave("outputs/figures/variant_classification_distribution.png", p2, width = 8, height = 6)

# --- 3. Mutations per patient (distribution histogram) ---
p3 <- mutation_counts_per_patient %>%
  ggplot(aes(x = mutation_count)) +
  geom_histogram(bins = 50, fill = "seagreen", color = "white") +
  scale_x_log10() +
  labs(title = "Distribution of Mutation Counts per Patient — TCGA-BRCA",
       x = "Mutation Count (log scale)", y = "Number of Patients") +
  theme_minimal()

ggsave("outputs/figures/mutations_per_patient.png", p3, width = 8, height = 6)

# --- 4. Mutation burden by clinical group (vital status) ---
patient_burden_clinical <- merged %>%
  distinct(patient_id, mutation_burden = NA, vital_status) %>%
  left_join(mutation_counts_per_patient, by = "patient_id") %>%
  filter(!is.na(mutation_count), !is.na(vital_status))

p4 <- patient_burden_clinical %>%
  ggplot(aes(x = vital_status, y = mutation_count, fill = vital_status)) +
  geom_boxplot() +
  scale_y_log10() +
  labs(title = "Mutation Burden by Vital Status — TCGA-BRCA",
       x = "Vital Status", y = "Mutation Count (log scale)") +
  theme_minimal() +
  theme(legend.position = "none")

ggsave("outputs/figures/mutation_burden_by_clinical_group.png", p4, width = 7, height = 6)

message("✅ Saved 4 figures to outputs/figures/")