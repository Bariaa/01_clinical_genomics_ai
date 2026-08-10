# Dataset Selection Report

## Objective
Document the rationale for assigning TCGA-BRCA, TCGA-LUAD, and TCGA-LAML to
primary (development), secondary (validation), and tertiary (stress test) roles
in evaluating **pipeline** — not model — generalizability.

## Rationale

### Primary dataset: TCGA-BRCA
Selected as the primary/development dataset due to its size (1088 patients) and
depth of clinical, histologic, and molecular annotation. Developing the pipeline
against the most feature-complete cohort first forces early handling of missing
data, subtype heterogeneity, and multiple biospecimen types, producing a more
robust pipeline before it is tested elsewhere.

### Secondary dataset: TCGA-LUAD
Selected as the secondary/validation dataset because it is structurally similar
to BRCA — both are solid tumors with comparable staging and histology frameworks —
while remaining an independent cohort from a different organ system. This tests
generalizability across organ systems without confounding the result with a
fundamentally different data schema.

### Tertiary dataset: TCGA-LAML
Selected as a stress test rather than a standard validation set. As a hematologic
malignancy, LAML lacks solid-tumor fields (staging, laterality, histologic subtype)
entirely. Applying the unmodified pipeline here tests whether it fails gracefully
(clear errors/flags) or silently (incorrect output without warning) when its
structural assumptions do not hold — a more rigorous test of true robustness than
same-structure validation alone.

## Explicit Non-Goals
- This design does **not** test statistical/predictive model generalizability.
- LAML results are not expected to "pass" in the same sense as LUAD; failure here
  is informative, not disqualifying, provided it is a controlled/logged failure.

## Decision Log

| Date | Decision | Rationale |
|---|---|---|
| 2026-08-10 | BRCA set as primary | Largest, most annotated cohort |
| 2026-08-10 | LUAD set as secondary | Structurally comparable solid tumor, independent cohort |
| 2026-08-10 | LAML set as stress test | Structurally distinct; tests failure-mode robustness |

## References
- Weinstein JN, et al. The Cancer Genome Atlas Pan-Cancer analysis project. Nat Genet. 2013.
- TRIPOD / STROBE reporting guidelines (adapted for pipeline-focused, not model-focused, reporting).