#!/usr/bin/env Rscript
# Run one comparison and feature type. Input files must already be prepared.
# Usage: Rscript 04_FinalModel.R COMPARISON FEATURE_TYPE INPUT_ROOT OUTPUT_ROOT
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
OUT_DIR <- file.path(RESULT_ROOT, "04_FinalModel")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
task_prefix <- if (comparison == "IBD_HC") "IBD" else "CDUC"
label_col <- if (comparison == "IBD_HC") "Group_IBD" else "Group_CDUC"
pos_class <- if (comparison == "IBD_HC") "IBD" else "CD"
neg_class <- if (comparison == "IBD_HC") "HC" else "UC"
TRAIN_DATA <- file.path(INPUT_DIR, paste0(task_prefix, "_train_data_use.txt"))

# Enforce the discovery cohort boundary.
check_cohorts <- function(cohorts) {
  allowed <- c("PRJNA398089", "PRJNA1265906", "PRJNA1045596", "PRJNA389280", "HRA007915")
  cohorts <- sub("^(ENA-|GSA-)", "", as.character(cohorts))
  stopifnot(length(cohorts) > 0, !anyNA(cohorts), all(cohorts %in% allowed))
}
library(caret)
library(randomForest)
library(ROCR)
set.seed(123)
OUT_DIR_MODEL <- OUT_DIR
customRF <- list(
  type    = "Classification",
  library = "randomForest",
  loop    = NULL
)
customRF$parameters <- data.frame(
  parameter = c("mtry", "ntree", "nodesize", "maxnodes"),
  class     = rep("numeric", 4),
  label     = c("mtry", "ntree", "nodesize", "maxnodes")
)
customRF$grid <- function(x, y, len = NULL, search = "grid") NULL
customRF$fit <- function(x, y, wts, param, lev, last, weights, classProbs, ...) {
  randomForest(
    x        = x,
    y        = y,
    mtry     = param$mtry,
    ntree    = param$ntree,
    nodesize = param$nodesize,
    maxnodes = param$maxnodes,
    ...
  )
}
customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL) {
  predict(modelFit, newdata)
}
customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL) {
  predict(modelFit, newdata, type = "prob")
}
customRF$sort <- function(x) x[order(x[, "mtry"]), ]
customRF$levels <- function(x) x$classes

get_control_first_levels <- function(y_factor, pos_class) {
  c(setdiff(levels(factor(y_factor)), pos_class)[1], pos_class)
}

get_auc <- function(prob_pos, y_factor, pos_class) {
  neg_class <- setdiff(levels(factor(y_factor)), pos_class)[1]
  pred_obj <- ROCR::prediction(
    predictions = prob_pos,
    labels = as.character(y_factor),
    label.ordering = c(neg_class, pos_class)
  )
  auc_obj <- ROCR::performance(pred_obj, "auc")
  as.numeric(auc_obj@y.values[[1]])
}

make_pos_summary <- function(pos_class) {
  force(pos_class)
  function(data, lev = NULL, model = NULL) {
    obs <- factor(data$obs, levels = lev)
    pred <- factor(data$pred, levels = lev)
    neg_class <- setdiff(lev, pos_class)[1]
    c(
      ROC = get_auc(data[[pos_class]], obs, pos_class),
      Sens = caret::sensitivity(pred, obs, positive = pos_class),
      Spec = caret::specificity(pred, obs, negative = neg_class)
    )
  }
}

train_final_rf <- function(data_path, label_col, pos_class, task_prefix) {

  dat <- read.table(
    data_path,
    header      = TRUE,
    sep         = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  check_cohorts(dat$Cohort)
  feat_cols <- setdiff(colnames(dat), c("Sample", "Cohort", label_col))

  y <- factor(dat[[label_col]])
  X <- dat[, feat_cols, drop = FALSE]

  y <- factor(y, levels = get_control_first_levels(y, pos_class))

  set.seed(2025)

  ctrl <- trainControl(
    method        = "repeatedcv",
    number        = 5,
    repeats       = 3,
    classProbs    = TRUE,
    summaryFunction = make_pos_summary(pos_class),
    savePredictions = "final"
  )

  p <- ncol(X)

  tunegrid <- expand.grid(
    .mtry     = 1:ceiling(sqrt(p)),
    .ntree    = c(301, 501, 801, 1001),
    .nodesize = c(20, 50, 100),
    .maxnodes = c(5, 10, 15, 20)
  )

  fit_cv <- caret::train(
    x         = X,
    y         = y,
    method    = customRF,
    metric    = "ROC",
    preProcess = c("center", "scale"),
    tuneGrid  = tunegrid,
    trControl = ctrl
  )

  best_par <- fit_cv$bestTune

  out_best <- file.path(
    OUT_DIR_MODEL,
    paste0(task_prefix, "_Final_RF_bestTune.csv")
  )

  write.csv(best_par, file = out_best, row.names = FALSE)

  out_tune <- file.path(
    OUT_DIR_MODEL,
    paste0(task_prefix, "_Final_RF_tuning_results.csv")
  )
  write.csv(fit_cv$results, file = out_tune, row.names = FALSE)

  set.seed(3030)

  final_rf <- randomForest(
    x        = X,
    y        = y,
    mtry     = best_par$mtry,
    ntree    = best_par$ntree,
    nodesize = best_par$nodesize,
    maxnodes = best_par$maxnodes,
    importance = TRUE
  )

  model_rds <- file.path(
    OUT_DIR_MODEL,
    paste0(task_prefix, "_Final_RF_model.rds")
  )
  saveRDS(final_rf, file = model_rds)

  imp_mat <- randomForest::importance(final_rf)

  imp_df <- data.frame(
    Feature = rownames(imp_mat),
    imp_mat,
    stringsAsFactors = FALSE
  )

  imp_df <- imp_df[order(imp_df[, "MeanDecreaseGini"], decreasing = TRUE), ]

  out_imp <- file.path(
    OUT_DIR_MODEL,
    paste0(task_prefix, "_Final_RF_feature_importance.csv")
  )
  write.csv(imp_df, file = out_imp, row.names = FALSE)

  prob_mat <- predict(final_rf, X, type = "prob")
  prob_df  <- data.frame(
    Sample   = dat$Sample,
    Cohort   = dat$Cohort,
    Label    = y,
    Prob_pos = prob_mat[, pos_class],
    stringsAsFactors = FALSE
  )

  out_prob <- file.path(
    OUT_DIR_MODEL,
    paste0(task_prefix, "_Final_RF_train_predictions.csv")
  )
  write.csv(prob_df, file = out_prob, row.names = FALSE)

  invisible(
    list(
      model   = final_rf,
      bestPar = best_par,
      imp     = imp_df,
      pred    = prob_df
    )
  )
}
train_final_rf(TRAIN_DATA, label_col, pos_class, task_prefix)
