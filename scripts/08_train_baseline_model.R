# 08_train_baseline_model.R
# Trains baseline models on the PRIMARY dataset (TCGA-BRCA) against multiple
# possible prediction targets:
#   - survival_status (alive/dead)
#   - risk_group (high/low, derived from mutation burden as a proxy — true
#     clinical risk_group is not available in GDC data)
#   - relapse_status (progression_or_recurrence)
#   - treatment_response
#   - long_vs_short_survival (median split of overall_survival_days, deceased only)
#
# IMPORTANT: patient_id is used only for tracking/matching — never a predictor.
#
# For each target, feasibility (sample size, class balance) is checked first.
# Targets with insufficient data are reported but NOT force-trained, since a
# model trained on too few events per class is not meaningful.

library(caret)
library(dplyr)
library(yaml)

config <- yaml::read_yaml("config/config.yaml")
MIN_EVENTS_PER_CLASS <- 20   # minimum cases needed in the smaller class to attempt training

ml_features <- read.csv("data/processed/ml_feature_table.csv", stringsAsFactors = FALSE)
clinical <- read.csv("data/processed/clinical_cleaned.csv", stringsAsFactors = FALSE)

patient_ids <- ml_features$patient_id
base_predictors <- c("age", "sex", "mutation_burden", "TP53_mutated",
                     "PIK3CA_mutated", "driver_gene_mutated")

# =========================================================================
# BUILD EACH CANDIDATE TARGET
# =========================================================================

targets <- list()

# --- 1. survival_status (alive/dead) ---
targets$survival_status <- ml_features %>%
  transmute(patient_id, target = factor(outcome, levels = c("alive", "dead")))

# --- 2. risk_group (PROXY: mutation burden median split — true risk_group unavailable) ---
med_burden <- median(ml_features$mutation_burden, na.rm = TRUE)
targets$risk_group <- ml_features %>%
  transmute(
    patient_id,
    target = factor(ifelse(mutation_burden >= med_burden, "high_risk", "low_risk"),
                    levels = c("low_risk", "high_risk"))
  )

# --- 3. relapse_status ---
relapse_raw <- clinical %>% select(patient_id, progression_or_recurrence)
targets$relapse_status <- relapse_raw %>%
  filter(!is.na(progression_or_recurrence), !progression_or_recurrence %in% c("not reported", "")) %>%
  transmute(patient_id, target = factor(progression_or_recurrence))

# --- 4. treatment_response ---
tr_raw <- clinical %>% select(patient_id, treatments_pharmaceutical_treatment_outcome)
targets$treatment_response <- tr_raw %>%
  filter(!is.na(treatments_pharmaceutical_treatment_outcome),
         !treatments_pharmaceutical_treatment_outcome %in% c("Not Reported", "")) %>%
  transmute(patient_id, target = factor(treatments_pharmaceutical_treatment_outcome))

# --- 5. long_vs_short_survival (deceased patients only, median split of survival days) ---
survival_days <- clinical %>%
  filter(vital_status == "dead", !is.na(overall_survival_days))
med_survival <- median(survival_days$overall_survival_days, na.rm = TRUE)
targets$long_vs_short_survival <- survival_days %>%
  transmute(
    patient_id,
    target = factor(ifelse(overall_survival_days >= med_survival, "long_survival", "short_survival"),
                    levels = c("short_survival", "long_survival"))
  )

# =========================================================================
# FEASIBILITY CHECK FOR EACH TARGET
# =========================================================================

feasibility <- lapply(names(targets), function(name) {
  t <- targets[[name]]
  class_counts <- table(t$target)
  min_class_n <- if (length(class_counts) >= 2) min(class_counts) else 0
  data.frame(
    target = name,
    n_available = nrow(t),
    n_classes = length(class_counts),
    min_class_n = min_class_n,
    feasible = length(class_counts) >= 2 && min_class_n >= MIN_EVENTS_PER_CLASS
  )
}) %>% bind_rows()

message("=== Target feasibility summary ===")
print(feasibility)

dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(feasibility, "outputs/tables/prediction_target_feasibility.csv", row.names = FALSE)

# =========================================================================
# TRAIN A MODEL FOR EACH FEASIBLE TARGET
# =========================================================================

dir.create(config$paths$outputs_models, recursive = TRUE, showWarnings = FALSE)
all_results <- list()

for (name in feasibility$target[feasibility$feasible]) {
  
  message("\n--- Training model for target: ", name, " ---")
  
  target_df <- targets[[name]]
  
  # Exclude mutation_burden specifically when the target is derived from it
  # (risk_group is a mutation_burden median-split proxy) — including it would
  # leak the label definition directly into the model as a predictor.
  predictors_for_this_target <- if (name == "risk_group") {
    setdiff(base_predictors, "mutation_burden")
  } else {
    base_predictors
  }
  
  model_data <- ml_features %>%
    select(patient_id, all_of(predictors_for_this_target)) %>%
    inner_join(target_df, by = "patient_id") %>%
    mutate(sex = as.factor(sex)) %>%
    select(-patient_id) %>%   # patient_id excluded from modeling — tracking only
    na.omit()
  set.seed(config$project$seed)
  
  model <- tryCatch({
    train(
      target ~ .,
      data = model_data,
      method = "glm",
      family = "binomial",
      metric = "ROC",
      trControl = trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = twoClassSummary)
    )
  }, error = function(e) {
    message("⚠️ Training failed for ", name, ": ", conditionMessage(e))
    NULL
  })
  
  if (!is.null(model)) {
    auc <- round(max(model$results$ROC), 4)
    message(name, ": AUC-ROC = ", auc, " (n=", nrow(model_data), ")")
    saveRDS(model, file.path(config$paths$outputs_models, paste0("model_", name, ".rds")))
    all_results[[name]] <- data.frame(target = name, n = nrow(model_data), auc_roc = auc)
  }
}

results_summary <- bind_rows(all_results)
write.csv(results_summary, "outputs/tables/model_results_by_target.csv", row.names = FALSE)

message("\n=== Final results across all feasible targets ===")
print(results_summary)
message("✅ Feasibility report: outputs/tables/prediction_target_feasibility.csv")
message("✅ Results summary: outputs/tables/model_results_by_target.csv")