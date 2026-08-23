# 08_train_baseline_model.R
#
# LIMITATION NOTE: This dataset does not have a strong native clinical
# endpoint. survival_status (vital_status) showed near-chance predictive
# power (AUC ~0.556, see prior analysis), relapse_status was entirely
# missing, and treatment_response was too sparse (9 records) to model.
#
# EXPLORATORY ENDPOINT: risk_group — a simple binary classification target
# created by median-splitting mutation_burden into "high_risk"/"low_risk".
# This is NOT a validated clinical risk classification; it is an exploratory
# proxy chosen specifically because no strong clinical endpoint exists in
# this cohort. mutation_burden itself is excluded from the feature set to
# avoid leaking the label's own definition into the model.
#
# patient_id is used only for tracking — never a model feature.
#
# Models: Logistic Regression, Random Forest, XGBoost (optional)
#
# Expected outputs:
#   outputs/tables/model_performance_summary.csv
#   outputs/tables/feature_importance.csv
#   outputs/figures/confusion_matrix.png
#   outputs/figures/roc_curve.png
#   outputs/models/baseline_model.rds   (R equivalent of .pkl — see note above)

library(caret)
library(dplyr)
library(pROC)
library(yaml)
library(ggplot2)

config <- yaml::read_yaml("config/config.yaml")
set.seed(config$project$seed)

dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/models", recursive = TRUE, showWarnings = FALSE)

# =========================================================================
# BUILD THE EXPLORATORY ENDPOINT
# =========================================================================
ml_features <- read.csv("data/processed/ml_feature_table.csv", stringsAsFactors = FALSE)

med_burden <- median(ml_features$mutation_burden, na.rm = TRUE)
ml_features <- ml_features %>%
  mutate(risk_group = factor(ifelse(mutation_burden >= med_burden, "high_risk", "low_risk"),
                             levels = c("low_risk", "high_risk")))

# Predictors: exclude patient_id (tracking only), outcome (unused here),
# mutation_burden (defines the label — would leak), and risk_group itself
predictor_cols <- c("age", "sex", "TP53_mutated", "PIK3CA_mutated", "driver_gene_mutated")

model_data <- ml_features %>%
  select(all_of(predictor_cols), risk_group) %>%
  mutate(sex = as.factor(sex)) %>%
  na.omit()

message("Exploratory endpoint: risk_group | n = ", nrow(model_data),
        " | predictors: ", paste(predictor_cols, collapse = ", "))

# =========================================================================
# TRAIN / TEST SPLIT
# =========================================================================
train_idx <- createDataPartition(model_data$risk_group, p = 0.8, list = FALSE)
train_data <- model_data[train_idx, ]
test_data <- model_data[-train_idx, ]

ctrl <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = twoClassSummary)

# =========================================================================
# TRAIN MODELS
# =========================================================================
models <- list()

models$logistic_regression <- train(
  risk_group ~ ., data = train_data, method = "glm", family = "binomial",
  metric = "ROC", trControl = ctrl
)

if (!requireNamespace("randomForest", quietly = TRUE)) install.packages("randomForest")
models$random_forest <- train(
  risk_group ~ ., data = train_data, method = "rf",
  metric = "ROC", trControl = ctrl, importance = TRUE
)

# XGBoost — optional, wrapped so failure doesn't break the pipeline
xgb_available <- requireNamespace("xgboost", quietly = TRUE) ||
  tryCatch({ install.packages("xgboost"); TRUE }, error = function(e) FALSE)

if (xgb_available) {
  models$xgboost <- tryCatch(
    train(risk_group ~ ., data = train_data, method = "xgbTree",
          metric = "ROC", trControl = ctrl, verbosity = 0),
    error = function(e) { message("⚠️ XGBoost training failed: ", conditionMessage(e)); NULL }
  )
  if (is.null(models$xgboost)) models$xgboost <- NULL
}

models <- Filter(Negate(is.null), models)

# =========================================================================
# EVALUATE EACH MODEL ON THE HELD-OUT TEST SET
# =========================================================================
performance_summary <- list()
roc_data_all <- list()

for (model_name in names(models)) {
  model <- models[[model_name]]
  
  pred_class <- predict(model, newdata = test_data)
  pred_prob <- predict(model, newdata = test_data, type = "prob")[, "high_risk"]
  
  cm <- confusionMatrix(pred_class, test_data$risk_group, positive = "high_risk")
  roc_obj <- roc(response = test_data$risk_group, predictor = pred_prob, levels = c("low_risk", "high_risk"), quiet = TRUE)
  
  performance_summary[[model_name]] <- data.frame(
    model = model_name,
    accuracy = round(cm$overall["Accuracy"], 4),
    precision = round(cm$byClass["Precision"], 4),
    recall = round(cm$byClass["Recall"], 4),
    f1_score = round(cm$byClass["F1"], 4),
    roc_auc = round(as.numeric(auc(roc_obj)), 4)
  )
  
  roc_data_all[[model_name]] <- data.frame(
    model = model_name,
    specificity = rev(roc_obj$specificities),
    sensitivity = rev(roc_obj$sensitivities)
  )
  
  # Save confusion matrix figure for the best-performing model later; store CM object for now
  assign(paste0("cm_", model_name), cm)
}

performance_df <- bind_rows(performance_summary)
write.csv(performance_df, "outputs/tables/model_performance_summary.csv", row.names = FALSE)
message("\n=== Model performance summary ===")
print(performance_df)

# =========================================================================
# FEATURE IMPORTANCE (from Random Forest, most standard for this)
# =========================================================================
if ("random_forest" %in% names(models)) {
  rf_importance <- varImp(models$random_forest)$importance
  feature_importance <- data.frame(
    feature = rownames(rf_importance),
    importance = rf_importance[, 1]
  ) %>% arrange(desc(importance))
  
  write.csv(feature_importance, "outputs/tables/feature_importance.csv", row.names = FALSE)
  message("\n=== Feature importance (Random Forest) ===")
  print(feature_importance)
}

# =========================================================================
# CONFUSION MATRIX FIGURE (best model by ROC-AUC)
# =========================================================================
best_model_name <- performance_df$model[which.max(performance_df$roc_auc)]
best_cm <- get(paste0("cm_", best_model_name))

cm_table <- as.data.frame(best_cm$table)

p_cm <- ggplot(cm_table, aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = Freq), size = 6) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(title = paste0("Confusion Matrix — ", best_model_name, " (risk_group)"),
       x = "Actual", y = "Predicted") +
  theme_minimal()

ggsave("outputs/figures/confusion_matrix.png", p_cm, width = 6, height = 5)

# =========================================================================
# ROC CURVE FIGURE (all models overlaid)
# =========================================================================
roc_combined <- bind_rows(roc_data_all)

p_roc <- ggplot(roc_combined, aes(x = 1 - specificity, y = sensitivity, color = model)) +
  geom_line(linewidth = 1) +
  geom_abline(linetype = "dashed", color = "grey50") +
  labs(title = "ROC Curves — risk_group (exploratory endpoint)",
       x = "1 - Specificity", y = "Sensitivity") +
  theme_minimal()

ggsave("outputs/figures/roc_curve.png", p_roc, width = 7, height = 6)

# =========================================================================
# SAVE BEST MODEL (R equivalent of .pkl — see note at top of script)
# =========================================================================
saveRDS(models[[best_model_name]], "outputs/models/baseline_model.rds")

message("\n✅ Best model: ", best_model_name,
        " (ROC-AUC = ", performance_df$roc_auc[performance_df$model == best_model_name], ")")
message("✅ All outputs saved: model_performance_summary.csv, feature_importance.csv, confusion_matrix.png, roc_curve.png, baseline_model.rds")