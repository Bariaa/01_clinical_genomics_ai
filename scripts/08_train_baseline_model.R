# 08_train_baseline_model.R
# Trains a baseline model on the PRIMARY dataset (TCGA-BRCA) feature table.
#
# IMPORTANT: patient_id is used only for tracking/matching results — it is
# explicitly excluded from the model's feature set and never used as a predictor.
#
# NOTE: A focused baseline predictor set is used (age, sex, mutation_burden,
# TP53_mutated, PIK3CA_mutated, driver_gene_mutated) rather than all 37 sparse
# gene/variant-class columns. With only ~150 "dead" events in 1097 patients,
# using all sparse binary gene flags caused non-convergence and a degenerate
# model that just predicted the majority class (Kappa ~ 0). This smaller,
# clinically motivated set is a more defensible starting baseline.

library(caret)
library(dplyr)
library(yaml)

config <- yaml::read_yaml("config/config.yaml")

ml_features <- read.csv("data/processed/ml_feature_table.csv", stringsAsFactors = FALSE)

# --- Explicitly separate identifier from predictors ---
patient_ids <- ml_features$patient_id   # kept aside, for tracking only

predictor_cols <- setdiff(names(ml_features), c("patient_id", "outcome"))

# Defensive check: guarantee patient_id can never leak into the model
stopifnot(!"patient_id" %in% predictor_cols)

message("Available predictor columns (patient_id excluded): ", length(predictor_cols))

# --- Prepare modeling data ---
model_data <- ml_features %>%
  select(all_of(predictor_cols), outcome) %>%
  filter(!is.na(age), !is.na(outcome)) %>%
  mutate(
    sex = as.factor(sex),
    outcome = factor(outcome, levels = c("alive", "dead"))
  )

# --- Reduce to a focused baseline feature set (avoids overfitting on sparse gene flags) ---
baseline_predictors <- c("age", "sex", "mutation_burden", "TP53_mutated",
                         "PIK3CA_mutated", "driver_gene_mutated")

model_data_final <- model_data %>%
  select(all_of(baseline_predictors), outcome) %>%
  na.omit()

message("Final baseline modeling dataset: ", nrow(model_data_final), " patients, ",
        length(baseline_predictors), " predictors: ", paste(baseline_predictors, collapse = ", "))

# --- Train baseline model (using AUC-ROC, more informative than accuracy for imbalanced classes) ---
set.seed(config$project$seed)

model <- train(
  outcome ~ .,
  data = model_data_final,
  method = "glm",
  family = "binomial",
  metric = "ROC",
  trControl = trainControl(
    method = "cv", number = 5,
    classProbs = TRUE,
    summaryFunction = twoClassSummary
  )
)

print(model)

message("AUC-ROC: ", round(max(model$results$ROC), 4),
        " (0.5 = no better than chance, 1.0 = perfect discrimination)")

# --- Save model + predictions with patient_id re-attached for tracking ---
dir.create(config$paths$outputs_models, recursive = TRUE, showWarnings = FALSE)
saveRDS(model, file.path(config$paths$outputs_models, "baseline_model.rds"))

predictions <- ml_features %>%
  filter(!is.na(age), !is.na(outcome)) %>%
  select(all_of(baseline_predictors), outcome, patient_id) %>%
  na.omit() %>%
  mutate(sex = as.factor(sex), outcome = factor(outcome, levels = c("alive", "dead")))

predictions$predicted <- predict(model, newdata = predictions)
predictions$predicted_prob_dead <- predict(model, newdata = predictions, type = "prob")[, "dead"]

write.csv(
  predictions %>% select(patient_id, outcome, predicted, predicted_prob_dead),
  "outputs/tables/baseline_model_predictions.csv",
  row.names = FALSE
)

message("✅ Baseline model trained on TCGA-BRCA (n=", nrow(model_data_final), ") and saved.")
message("✅ Predictions (with patient_id for tracking only) saved to outputs/tables/baseline_model_predictions.csv")