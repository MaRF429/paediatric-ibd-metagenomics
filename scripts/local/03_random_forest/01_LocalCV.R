#!/usr/bin/env Rscript
# Run one comparison and feature type. Input files must already be prepared.
# Usage: Rscript 01_LocalCV.R COMPARISON FEATURE_TYPE INPUT_ROOT OUTPUT_ROOT
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
OUT_DIR <- file.path(RESULT_ROOT, "01_LocalCV")
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

customRF$grid <- function(x, y, len = NULL, search = "grid") {

  NULL
}

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

tune_model <- function(data_ml, y_factor, pos_class, num_repeat) {

  p <- ncol(data_ml)

  set.seed(2021 + num_repeat)
  y_factor <- factor(y_factor, levels = get_control_first_levels(y_factor, pos_class))

  ctrl <- trainControl(
    method        = "repeatedcv",
    number        = 5,
    repeats       = 1,
    classProbs    = TRUE,
    summaryFunction = make_pos_summary(pos_class),
    savePredictions = "final"
  )

  tunegrid <- expand.grid(
    .mtry     = 1:ceiling(sqrt(p)),
    .ntree    = c(301, 501, 801, 1001),
    .nodesize = c(50, 100, 150),
    .maxnodes = c(5, 10, 15, 20)
  )

  y_factor <- factor(y_factor, levels = get_control_first_levels(y_factor, pos_class))

  rf_fit <- caret::train(
    x         = data_ml,
    y         = y_factor,
    method    = customRF,
    metric    = "ROC",
    preProcess = c("center", "scale"),
    tuneGrid  = tunegrid,
    trControl = ctrl
  )

  best_rf <- rf_fit$bestTune

  return(best_rf)
}

get_self_cv <- function(data_s, label_col, best_rf, pos_class, num_repeat) {

  set.seed(2021 + num_repeat)

  y_all <- factor(data_s[[label_col]])

  y_all <- factor(y_all, levels = get_control_first_levels(y_all, pos_class))

  X_all <- data_s[, setdiff(colnames(data_s), label_col), drop = FALSE]

  folds <- createFolds(y = y_all, k = 5, list = TRUE, returnTrain = FALSE)

  metrics <- matrix(NA, nrow = 5, ncol = 12)
  colnames(metrics) <- c(
    "Sensitivity", "Specificity", "Pos Pred Value", "Neg Pred Value",
    "Precision", "Recall", "F1", "Prevalence", "Detection Rate",
    "Detection Prevalence", "Balanced Accuracy", "AUC"
  )

  for (j in seq_along(folds)) {
    idx_te  <- folds[[j]]
    idx_tr  <- setdiff(seq_len(nrow(data_s)), idx_te)

    X_tr <- X_all[idx_tr, , drop = FALSE]
    y_tr <- y_all[idx_tr]

    X_te <- X_all[idx_te, , drop = FALSE]
    y_te <- y_all[idx_te]

    set.seed(518 + j)

    rf_fit <- randomForest(
      x        = X_tr,
      y        = y_tr,
      mtry     = best_rf$mtry,
      ntree    = best_rf$ntree,
      nodesize = best_rf$nodesize,
      maxnodes = best_rf$maxnodes,
      importance = TRUE
    )

    pred_class <- predict(rf_fit, X_te)

    cm <- confusionMatrix(
      data      = factor(pred_class, levels = levels(y_all)),
      reference = y_te,
      mode      = "everything",
      positive  = pos_class
    )

    metrics[j, 1:11] <- as.numeric(cm$byClass)

    pred_prob <- predict(rf_fit, X_te, type = "prob")

    prob_pos  <- pred_prob[, pos_class]

    metrics[j, 12] <- get_auc(prob_pos, y_te, pos_class)
  }

  best_index <- which(metrics[, "AUC"] == max(metrics[, "AUC"], na.rm = TRUE))[1]

  best_model <- list(
    best_fold  = best_index,
    metrics    = metrics
  )

  return(best_model)
}

  dat <- read.table(
    TRAIN_DATA,
    header      = TRUE,
    sep         = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  check_cohorts(dat$Cohort)


  feat_cols <- setdiff(colnames(dat), c("Sample", "Cohort", label_col))

  cohorts <- sort(unique(dat$Cohort))

  for (s in cohorts) {

    dat_s_all <- dat[dat$Cohort == s, , drop = FALSE]

    data_s <- dat_s_all[, c(feat_cols, label_col), drop = FALSE]

    Group_s <- factor(data_s[[label_col]])

    if (length(unique(Group_s)) < 2) {

      next
    }

    model_s      <- matrix(NA, nrow = 20, ncol = 4)
    colnames(model_s) <- c("mtry", "ntree", "nodesize", "maxnodes")

    s_matrix     <- matrix(NA, nrow = 20 * 5, ncol = 12)
    colnames(s_matrix) <- c(
      "Sensitivity","Specificity","Pos Pred Value","Neg Pred Value","Precision",
      "Recall","F1","Prevalence","Detection Rate","Detection Prevalence",
      "Balanced Accuracy","AUC"
    )

    s_matrix_max <- matrix(NA, nrow = 20, ncol = 12)
    colnames(s_matrix_max) <- colnames(s_matrix)

    for (i in 1:20) {

      data_ml <- data_s[, feat_cols, drop = FALSE]
      best_rf <- tune_model(
        data_ml  = data_ml,
        y_factor = Group_s,
        pos_class = pos_class,
        num_repeat = i
      )

      model_s[i, ] <- c(best_rf$mtry, best_rf$ntree, best_rf$nodesize, best_rf$maxnodes)

      s_AUC <- get_self_cv(
        data_s    = data_s,
        label_col = label_col,
        best_rf   = best_rf,
        pos_class = pos_class,
        num_repeat = i
      )

      idx_start <- (i - 1) * 5 + 1
      idx_end   <- i * 5
      s_matrix[idx_start:idx_end, ] <- s_AUC$metrics

      s_matrix_max[i, ] <- s_AUC$metrics[s_AUC$best_fold, ]
    }

    mean_auc_all <- mean(s_matrix[, "AUC"], na.rm = TRUE)
    mean_auc_max <- mean(s_matrix_max[, "AUC"], na.rm = TRUE)

    prefix <- paste0(task_prefix, "_LocalCV_RF_", s)

    write.csv(
      s_matrix,
      file = file.path(OUT_DIR, paste0(prefix, "_selfcv_repeat_metric.csv")),
      row.names = FALSE
    )

    write.csv(
      s_matrix_max,
      file = file.path(OUT_DIR, paste0(prefix, "_selfcv_repeat_max.csv")),
      row.names = FALSE
    )

    write.csv(
      model_s,
      file = file.path(OUT_DIR, paste0(prefix, "_selfcv_repeat_bestTune.csv")),
      row.names = FALSE
    )
  }
