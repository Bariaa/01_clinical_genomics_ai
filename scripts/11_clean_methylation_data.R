# 11_clean_methylation_data.R
# Clean and QC DNA methylation beta values for the PRIMARY dataset (TCGA-BRCA).
#
# Steps:
#   1. Load the reduced beta matrix from script 10
#   2. Assess per-probe and per-sample missingness
#   3. Filter probes exceeding a missingness threshold
#   4. Filter samples exceeding a missingness threshold
#   5. Impute remaining NAs (probe-wise median)
#   6. Summarise probe characteristics (mean beta, SD, methylation class)
#   7. Save cleaned matrix and summary tables
#
# Outputs:
#   data/processed/methylation_cleaned.csv
#   outputs/tables/methylation_probe_summary.csv
#   outputs/tables/methylation_missingness_summary.csv

library(data.table)
library(yaml)

config <- yaml::read_yaml("config/config.yaml")
primary_id <- config$datasets$primary$name

MAX_PROBE_MISSING  <- 0.20   # drop probes missing in >20% of samples
MAX_SAMPLE_MISSING <- 0.20   # drop samples missing >20% of probes

beta <- readRDS(file.path(config$paths$data_raw, primary_id,
                          paste0(primary_id, "_methylation_betas.rds")))
message("Loaded beta matrix: ", paste(dim(beta), collapse = " x "))

# --- 2. Missingness assessment (before filtering) ---
probe_missing_before  <- rowMeans(is.na(beta))
sample_missing_before <- colMeans(is.na(beta))

missingness_summary <- data.frame(
  stage = "before_filtering",
  n_probes = nrow(beta),
  n_samples = ncol(beta),
  overall_pct_missing = round(100 * mean(is.na(beta)), 3),
  probes_over_threshold = sum(probe_missing_before > MAX_PROBE_MISSING),
  samples_over_threshold = sum(sample_missing_before > MAX_SAMPLE_MISSING),
  median_probe_pct_missing = round(100 * median(probe_missing_before), 3),
  median_sample_pct_missing = round(100 * median(sample_missing_before), 3)
)

message("Probes over missingness threshold: ", missingness_summary$probes_over_threshold)
message("Samples over missingness threshold: ", missingness_summary$samples_over_threshold)

# --- 3 & 4. Filter ---
keep_probes  <- probe_missing_before  <= MAX_PROBE_MISSING
keep_samples <- sample_missing_before <= MAX_SAMPLE_MISSING

beta_filt <- beta[keep_probes, keep_samples, drop = FALSE]
message("After filtering: ", paste(dim(beta_filt), collapse = " x "))

# --- 5. Impute remaining NAs with probe-wise median ---
n_na_before_impute <- sum(is.na(beta_filt))
if (n_na_before_impute > 0) {
  message("Imputing ", n_na_before_impute, " remaining NA values (probe-wise median)...")
  probe_medians <- apply(beta_filt, 1, median, na.rm = TRUE)
  na_idx <- which(is.na(beta_filt), arr.ind = TRUE)
  beta_filt[na_idx] <- probe_medians[na_idx[, "row"]]
}
message("NAs remaining after imputation: ", sum(is.na(beta_filt)))

# --- Missingness summary (after) ---
missingness_summary <- rbind(
  missingness_summary,
  data.frame(
    stage = "after_filtering_and_imputation",
    n_probes = nrow(beta_filt),
    n_samples = ncol(beta_filt),
    overall_pct_missing = round(100 * mean(is.na(beta_filt)), 3),
    probes_over_threshold = 0L,
    samples_over_threshold = 0L,
    median_probe_pct_missing = 0,
    median_sample_pct_missing = 0
  )
)

# --- 6. Probe summary / annotation ---
# Methylation class follows the conventional beta-value cut-offs:
#   < 0.3 = hypomethylated, 0.3-0.7 = intermediate, > 0.7 = hypermethylated
probe_mean <- rowMeans(beta_filt)
probe_sd   <- apply(beta_filt, 1, sd)

probe_summary <- data.frame(
  probe_id = rownames(beta_filt),
  mean_beta = round(probe_mean, 4),
  sd_beta = round(probe_sd, 4),
  min_beta = round(apply(beta_filt, 1, min), 4),
  max_beta = round(apply(beta_filt, 1, max), 4),
  pct_missing_before_filter = round(100 * probe_missing_before[rownames(beta_filt)], 3),
  methylation_class = cut(probe_mean,
                          breaks = c(-Inf, 0.3, 0.7, Inf),
                          labels = c("hypomethylated", "intermediate", "hypermethylated")),
  stringsAsFactors = FALSE
)
probe_summary <- probe_summary[order(-probe_summary$sd_beta), ]

message("\nMethylation class distribution:")
print(table(probe_summary$methylation_class))

# --- 7. Save ---
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

# Cleaned matrix: probes as rows, samples as columns, probe_id as first column
out <- data.table::data.table(probe_id = rownames(beta_filt),
                              as.data.frame(beta_filt))
data.table::fwrite(out, "data/processed/methylation_cleaned.csv")

write.csv(probe_summary, "outputs/tables/methylation_probe_summary.csv", row.names = FALSE)
write.csv(missingness_summary, "outputs/tables/methylation_missingness_summary.csv", row.names = FALSE)

# Keep an .rds alongside for fast downstream loading (gitignored)
saveRDS(beta_filt, file.path(config$paths$data_processed,
                             paste0(primary_id, "_methylation_cleaned.rds")))

message("\n✅ Saved: data/processed/methylation_cleaned.csv (",
        round(file.size("data/processed/methylation_cleaned.csv")/1024^2, 1), " MB)")
message("✅ Saved: outputs/tables/methylation_probe_summary.csv")
message("✅ Saved: outputs/tables/methylation_missingness_summary.csv")
print(missingness_summary)