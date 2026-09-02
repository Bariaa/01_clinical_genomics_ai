# LUAD Generalizability Test Report

## Objective
Test whether the pipeline built for TCGA-BRCA (primary dataset) runs
unmodified on TCGA-LUAD (secondary dataset), and identify what does/does
not generalize.

## Can the scripts run on the second dataset?
**Yes.** Scripts 02 (clinical cleaning), 03 (mutation processing), and 04
(ID matching) all ran successfully on TCGA-LUAD without any code changes,
using the same column-name assumptions, ID-parsing logic, and join keys
built for BRCA.

## Which parts work without changes?
- **Patient ID extraction** (first 12 characters of TCGA barcode) — identical
  format across both cancer types, worked unmodified.
- **Clinical cleaning logic** (janitor::clean_names, list-column flattening,
  deduplication) — ran cleanly, same GDC schema conventions apply.
- **Mutation processing** (Hugo_Symbol, Chromosome, Variant_Classification,
  Variant_Type, Tumor_Sample_Barcode fields) — same standard MAF schema
  across both cancer types, no adjustment needed.
- **ID matching/merge logic** — worked at an even higher match rate (95.21%
  vs 88.16% for BRCA).
- **survival_status field** (vital_status) — present and usable, same as BRCA.

## Which parts need adjustment?
- **Gene-based features** (TP53_mutated, PIK3CA_mutated, driver_gene list)
  were curated specifically around BRCA biology. In LUAD, PIK3CA does not
  appear in the top 10 mutated genes at all, and TP53 ranks 7th rather than
  1st. A LUAD-specific driver gene list (e.g., KRAS, EGFR, STK11, KEAP1 —
  well-established LUAD drivers) would be needed for a meaningful
  LUAD-specific model; the BRCA-curated gene features would carry much
  weaker signal here.
- **risk_group proxy threshold** (mutation_burden median split) would need
  recalculating for LUAD specifically, since its overall mutation burden
  distribution is substantially higher than BRCA's.

## Are the column names different?
Mostly no — LUAD's cleaned clinical table has 101 columns vs BRCA's 103,
a small difference (2 columns), suggesting near-identical GDC clinical
schemas between the two cancer types. No renaming was required for the
pipeline to run.

## Are the patient/sample ID formats different?
No. Both use the standard TCGA barcode format
(`TCGA-XX-XXXX-...`), with the first 12 characters identifying the patient
consistently across cancer types. The same substring-based ID extraction
logic worked unmodified.

## Are the available clinical endpoints different?
Partially. `vital_status` is present and usable in both. Critically, LUAD
has a substantially more balanced outcome distribution — 32.1% deceased
(188/585) vs BRCA's 13.86% — which likely makes `survival_status` a *more*
learnable classification target in LUAD than it was in BRCA, where severe
class imbalance was a contributing factor to the near-chance baseline
result. `relapse_status` and `treatment_response` were not re-tested here
but are drawn from the same GDC fields and would need the same feasibility
check performed for BRCA, since sparsity is likely dataset-specific.

## What should be improved to make the workflow more general?
1. **Externalize gene lists to config**: driver genes and genes-of-interest
   should be defined per-dataset in `config.yaml` rather than hardcoded in
   scripts, since disease-relevant genes differ substantially by cancer type.
2. **Recompute proxy thresholds per dataset**: any derived label (like
   risk_group's median split) should be calculated fresh for each dataset,
   not inherited from the primary dataset's distribution.
3. **Automate the feasibility check** (already built for BRCA's target
   selection in script 08) as a standard step for any new dataset, before
   attempting to train models on clinical endpoints.
4. **Parameterize column-count/schema comparison** as an automatic QC step
   (e.g., diff BRCA vs LUAD column names) rather than manual inspection,
   to catch schema drift early for future datasets.

