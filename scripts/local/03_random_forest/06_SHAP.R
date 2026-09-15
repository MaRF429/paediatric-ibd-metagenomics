#!/usr/bin/env Rscript
# Run one comparison and feature type. Input files must already be prepared.
# Usage: Rscript 06_SHAP.R COMPARISON FEATURE_TYPE INPUT_ROOT OUTPUT_ROOT
# Inputs and outputs are isolated under COMPARISON/FEATURE_TYPE.
options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("Expected comparison, feature type, input root and output root.")
comparison <- args[[1]]
feature_type <- args[[2]]
stopifnot(comparison %in% c("IBD_HC", "CD_UC"),
          feature_type %in% c("species", "KO", "EggNOG", "pathway", "combined"))

INPUT_DIR <- file.path(args[[3]], comparison, feature_type)
RESULT_ROOT <- file.path(args[[4]], comparison, feature_type)
OUT_DIR <- file.path(RESULT_ROOT, "06_SHAP")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
task_prefix <- if (comparison == "IBD_HC") "IBD" else "CDUC"
label_col <- if (comparison == "IBD_HC") "Group_IBD" else "Group_CDUC"
pos_class <- if (comparison == "IBD_HC") "IBD" else "CD"
TRAIN_DATA <- file.path(INPUT_DIR, paste0(task_prefix, "_train_data_use.txt"))
MODEL_PATH <- file.path(RESULT_ROOT, "04_FinalModel", paste0(task_prefix, "_Final_RF_model.rds"))

# Enforce the discovery cohort boundary.
check_cohorts <- function(cohorts) {
  allowed <- c("PRJNA398089", "PRJNA1265906", "PRJNA1045596", "PRJNA389280", "HRA007915")
  cohorts <- sub("^(ENA-|GSA-)", "", as.character(cohorts))
  stopifnot(length(cohorts) > 0, !anyNA(cohorts), all(cohorts %in% allowed))
}

library(randomForest)
library(fastshap)
set.seed(123)
dat <- read.table(TRAIN_DATA, header = TRUE, sep = "\t", check.names = FALSE)
check_cohorts(dat$Cohort)
feat_cols <- setdiff(colnames(dat), c("Sample", "Cohort", label_col))
X <- dat[, feat_cols, drop = FALSE]
rf_model <- readRDS(MODEL_PATH)
model_features <- row.names(rf_model$importance)
X <- X[, model_features, drop = FALSE]
pred_fun <- function(object, newdata) {
  predict(object, newdata, type = "prob")[, pos_class]
}
shap_values <- fastshap::explain(rf_model, X = X,
                                pred_wrapper = pred_fun, nsim = 500)
write.csv(shap_values, file.path(OUT_DIR, "SHAP_values.csv"))
importance <- data.frame(Feature = colnames(shap_values),
                          MeanAbsSHAP = colMeans(abs(shap_values)))
importance <- importance[order(importance$MeanAbsSHAP, decreasing = TRUE), ]
importance$Rank <- seq_len(nrow(importance))
write.csv(importance, file.path(OUT_DIR, "SHAP_importance_table.csv"), row.names = FALSE)
