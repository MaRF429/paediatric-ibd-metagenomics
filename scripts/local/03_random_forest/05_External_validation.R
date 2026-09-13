#!/usr/bin/env Rscript
# Run one comparison and feature type. Input files must already be prepared.
# Usage: Rscript 05_External_validation.R COMPARISON FEATURE_TYPE INPUT_ROOT OUTPUT_ROOT
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
OUT_DIR <- file.path(RESULT_ROOT, "05_External_validation")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
task_prefix <- if (comparison == "IBD_HC") "IBD" else "CDUC"
label_col <- if (comparison == "IBD_HC") "Group_IBD" else "Group_CDUC"
pos_class <- if (comparison == "IBD_HC") "IBD" else "CD"
neg_class <- if (comparison == "IBD_HC") "HC" else "UC"
TEST_DATA <- file.path(INPUT_DIR, paste0(task_prefix, "_test_data_use.txt"))
MODEL_PATH <- file.path(RESULT_ROOT, "04_FinalModel", paste0(task_prefix, "_Final_RF_model.rds"))

# Enforce the independent-validation cohort boundary.
check_cohorts <- function(cohorts) {
  allowed <- c("PRJEB76677", "SRP057027", "PRJNA759642", "PRJNA922068")
  cohorts <- sub("^(ENA-|GSA-)", "", as.character(cohorts))
  stopifnot(length(cohorts) > 0, !anyNA(cohorts), all(cohorts %in% allowed))
}
library(caret)
library(randomForest)
library(ROCR)
set.seed(123)
get_auc <- function(prob_pos, y_factor, pos_class) {

  neg_class <- setdiff(levels(factor(y_factor)), pos_class)[1]
  pred_obj <- ROCR::prediction(
    predictions = prob_pos,
    labels = as.character(y_factor),
    label.ordering = c(neg_class, pos_class)
  )
  auc_obj  <- ROCR::performance(pred_obj, "auc")
  return(as.numeric(auc_obj@y.values[[1]]))
}
get_control_first_levels <- function(y_factor, pos_class) {
  c(setdiff(levels(factor(y_factor)), pos_class)[1], pos_class)
}
metric_names <- c(
  "Sensitivity","Specificity","Pos Pred Value","Neg Pred Value",
  "Precision","Recall","F1","Prevalence","Detection Rate",
  "Detection Prevalence","Balanced Accuracy","AUC"
)

external_test_one_cohort <- function(final_rf, dat_cohort, label_col, pos_class) {

  feat_cols <- setdiff(colnames(dat_cohort), c("Sample","Cohort", label_col))
  X <- dat_cohort[, feat_cols, drop = FALSE]
  y <- factor(dat_cohort[[label_col]])
  y <- factor(y, levels = get_control_first_levels(y, pos_class))

  if (length(unique(y)) < 2) {
    warning("Required input is unavailable or invalid.")
    return(rep(NA, length(metric_names)))
  }

  pred_class <- predict(final_rf, X)
  pred_prob  <- predict(final_rf, X, type = "prob")[, pos_class]

  cm <- confusionMatrix(
    data      = pred_class,
    reference = y,
    mode      = "everything",
    positive  = pos_class
  )

  m_vec <- numeric(length(metric_names))
  m_vec[1:11] <- as.numeric(cm$byClass)
  m_vec[12]   <- get_auc(pred_prob, y, pos_class)

  m_vec
}

# Locked-model evaluation only: no training or tuning on validation data.
dat_te <- read.table(TEST_DATA, header = TRUE, sep = "\t",
                     check.names = FALSE, stringsAsFactors = FALSE)
check_cohorts(dat_te$Cohort)
final_rf <- readRDS(MODEL_PATH)
val_cohorts <- sort(unique(dat_te$Cohort))
ext_metrics_mat <- matrix(NA, nrow = length(val_cohorts) + 1,
                          ncol = length(metric_names),
                          dimnames = list(c(val_cohorts, "All_validation"), metric_names))
for (co in val_cohorts) {
  ext_metrics_mat[co, ] <- external_test_one_cohort(
    final_rf, dat_te[dat_te$Cohort == co, , drop = FALSE], label_col, pos_class)
}
ext_metrics_mat["All_validation", ] <- external_test_one_cohort(
  final_rf, dat_te, label_col, pos_class)
write.csv(ext_metrics_mat, file.path(OUT_DIR, paste0(task_prefix, "_External_FinalModel_metrics.csv")))
