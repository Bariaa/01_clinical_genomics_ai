# 15_build_epigenomic_feature_table.R
# Build the AI-ready multi-omics feature table.

# 15_build_epigenomic_feature_table.R
# Build the AI-ready multi-omics feature table.
#
# The integrated dataset holds 17,095 methylation features for 649 patients
# (~26 features per patient). Modelling on that directly would overfit, so
# the methylation layer is reduced to three complementary representations:
#   1. Top-variance CpG probes (individual, interpretable)
#   2. Principal components (global structure, compact)
#   3. Summary statistics (mean beta, hypo/hyper probe counts)
#
# LEAKAGE NOTE: all three reductions are UNSUPERVISED — variance ranking and
# PCA never see the outcome variable. This matters: in Month 1 a supervised
# feature (mutation_burden) was accidentally used to both define a label and
# predict it, producing a spurious AUC of 1.0.
#
# Outputs:
#   data/processed/multiomics_feature_table.csv
#   outputs/tables/epigenomic_feature_summary.csv

library(dplyr)
library(yaml)

config <- yaml::read_yaml("config/config.yaml")
proc_dir <- config$paths$data_processed
set.seed(config$project$seed)

N_TOP_PROBES <- 500
N_PCS <- 20

integrated <- readRDS(file.path(proc_dir, "multiomics_integrated.rds"))
meth_cols <- grep("^cg", names(integrated), value = TRUE)
beta_mat <- as.matrix(integrated[, meth_cols])
rownames(beta_mat) <- integrated$patient_id

message("Input: ", nrow(beta_mat), " patients x ", ncol(beta_mat), " methylation features")

# --- 1. Top-variance probes (unsupervised) ---
probe_sd <- apply(beta_mat, 2, sd)
top_probes <- names(sort(probe_sd, decreasing = TRUE))[seq_len(N_TOP_PROBES)]
probe_features <- as.data.frame(beta_mat[, top_probes, drop = FALSE])
names(probe_features) <- paste0("meth_", names(probe_features))
message("Selected top ", N_TOP_PROBES, " variable probes")

# --- 2. Principal components (unsupervised) ---
message("Computing ", N_PCS, " principal components...")
pca <- prcomp(beta_mat, center = TRUE, scale. = FALSE)
pc_features <- as.data.frame(pca$x[, seq_len(N_PCS), drop = FALSE])
names(pc_features) <- paste0("meth_PC", seq_len(N_PCS))
pc_var <- round(100 * summary(pca)$importance[2, seq_len(N_PCS)], 2)
message("Variance captured by ", N_PCS, " PCs: ",
        round(sum(summary(pca)$importance[2, seq_len(N_PCS)]) * 100, 1), "%")

# --- 3. Summary statistics ---
summary_features <- data.frame(
  meth_mean_beta        = rowMeans(beta_mat),
  meth_sd_beta          = apply(beta_mat, 1, sd),
  meth_median_beta      = apply(beta_mat, 1, median),
  meth_n_hypomethylated = rowSums(beta_mat < 0.3),
  meth_n_intermediate   = rowSums(beta_mat >= 0.3 & beta_mat <= 0.7),
  meth_n_hypermethylated= rowSums(beta_mat > 0.7)
)
message("Computed ", ncol(summary_features), " summary features")

# --- 4. Clinical and genomic features from Month 1 ---
clinical_genomic_cols <- setdiff(names(integrated), c(meth_cols, "patient_id"))
cg_features <- integrated[, clinical_genomic_cols, drop = FALSE]
message("Clinical + genomic features carried forward: ", ncol(cg_features))
message("  ", paste(clinical_genomic_cols, collapse = ", "))

# --- 5. Assemble the feature table ---
feature_table <- bind_cols(
  data.frame(patient_id = integrated$patient_id, stringsAsFactors = FALSE),
  cg_features,
  summary_features,
  pc_features,
  probe_features
)

message("\nFeature table: ", nrow(feature_table), " patients x ",
        ncol(feature_table), " columns")

# --- 6. Feature inventory (documents what went in and why) ---
feature_summary <- bind_rows(
  data.frame(
    feature_group = "clinical",
    n_features = sum(clinical_genomic_cols %in%
                       c("age_at_index", "sex_at_birth", "primary_diagnosis",
                         "vital_status", "overall_survival_days", "survival_event")),
    selection_method = "carried forward from Month 1",
    supervised = "no",
    notes = "patient-level clinical annotation"
  ),
  data.frame(
    feature_group = "genomic",
    n_features = sum(clinical_genomic_cols %in%
                       c("mutation_burden", "TP53_mutated", "PIK3CA_mutated")),
    selection_method = "carried forward from Month 1",
    supervised = "no",
    notes = "mutation burden and driver gene status"
  ),
  data.frame(
    feature_group = "methylation_summary",
    n_features = ncol(summary_features),
    selection_method = "computed across all 17,095 probes",
    supervised = "no",
    notes = "mean/median/sd beta, hypo/intermediate/hyper probe counts"
  ),
  data.frame(
    feature_group = "methylation_pcs",
    n_features = ncol(pc_features),
    selection_method = paste0("PCA, first ", N_PCS, " components"),
    supervised = "no",
    notes = paste0("captures ", round(sum(summary(pca)$importance[2, seq_len(N_PCS)]) * 100, 1),
                   "% of methylation variance")
  ),
  data.frame(
    feature_group = "methylation_probes",
    n_features = ncol(probe_features),
    selection_method = paste0("top ", N_TOP_PROBES, " by standard deviation"),
    supervised = "no",
    notes = "individual CpG sites, biologically interpretable"
  )
)
feature_summary$pct_of_total <- round(100 * feature_summary$n_features /
                                        sum(feature_summary$n_features), 2)

# --- 7. Save ---
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(feature_table, "data/processed/multiomics_feature_table.csv", row.names = FALSE)
write.csv(feature_summary, "outputs/tables/epigenomic_feature_summary.csv", row.names = FALSE)
saveRDS(feature_table, file.path(proc_dir, "multiomics_feature_table.rds"))

# PC variance table — useful for the report
pc_variance <- data.frame(
  component = names(pc_var),
  pct_variance = as.numeric(pc_var),
  cumulative_pct = round(cumsum(as.numeric(pc_var)), 2)
)
write.csv(pc_variance, "outputs/tables/methylation_pc_variance.csv", row.names = FALSE)

message("\n✅ Saved: data/processed/multiomics_feature_table.csv (",
        round(file.size("data/processed/multiomics_feature_table.csv")/1024^2, 1), " MB)")
message("✅ Saved: outputs/tables/epigenomic_feature_summary.csv")
message("✅ Saved: outputs/tables/methylation_pc_variance.csv")
print(feature_summary)

