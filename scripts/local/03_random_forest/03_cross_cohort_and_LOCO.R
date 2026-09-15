#!/usr/bin/env Rscript
# Run one comparison and feature type. Input files must already be prepared.
# Usage: Rscript 03_cross_cohort_and_LOCO.R COMPARISON FEATURE_TYPE INPUT_ROOT OUTPUT_ROOT
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
OUT_DIR <- file.path(RESULT_ROOT, "03_cross_cohort_and_LOCO")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
task_prefix <- if (comparison == "IBD_HC") "IBD" else "CDUC"
label_col <- if (comparison == "IBD_HC") "Group_IBD" else "Group_CDUC"
pos_class <- if (comparison == "IBD_HC") "IBD" else "CD"
neg_class <- if (comparison == "IBD_HC") "HC" else "UC"
TRAIN_DATA <- file.path(INPUT_DIR, paste0(task_prefix, "_train_data_use.txt"))
LOCALCV_DIR <- file.path(RESULT_ROOT, "01_LocalCV")

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


OUT_DIR_CC <- file.path(OUT_DIR, "CrossCohort_RF")
OUT_DIR_LOCO <- file.path(OUT_DIR, "LOCO_RF")
dir.create(OUT_DIR_CC, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR_LOCO, recursive = TRUE, showWarnings = FALSE)
get_representative_params <- function(task_prefix, cohort) {
  file_para <- file.path(
    LOCALCV_DIR,
    paste0(task_prefix, "_LocalCV_RF_", cohort, "_selfcv_repeat_bestTune.csv")
  )
  para_df <- read.csv(file_para, header = TRUE, check.names = FALSE)
  if (nrow(para_df) < 1) {
    stop("Required input is unavailable or invalid.", file_para)
  }
  rep_row <- para_df[1, ]
  list(
    mtry     = rep_row[["mtry"]],
    ntree    = rep_row[["ntree"]],
    nodesize = rep_row[["nodesize"]],
    maxnodes = rep_row[["maxnodes"]]
  )
}

get_auc <- function(prob_pos, y_factor, pos_class) {

  neg_class <- setdiff(levels(y_factor), pos_class)[1]
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

metric_names <- c(
  "Sensitivity","Specificity","Pos Pred Value","Neg Pred Value",
  "Precision","Recall","F1","Prevalence","Detection Rate",
  "Detection Prevalence","Balanced Accuracy","AUC"
)

run_cross_cohort <- function(dat, task_prefix, label_col, pos_class, out_dir) {

  feat_cols <- setdiff(colnames(dat), c("Sample","Cohort", label_col))
  cohorts   <- sort(unique(dat$Cohort))

  all_pairs <- character()
  for (tr in cohorts) {
    for (te in cohorts) {
      if (te != tr) {
        all_pairs <- c(all_pairs, paste(tr, te, sep = "_"))
      }
    }
  }
  metrics_mat <- matrix(NA, nrow = length(all_pairs), ncol = length(metric_names))
  rownames(metrics_mat) <- all_pairs
  colnames(metrics_mat) <- metric_names

  for (tr in cohorts) {

    dat_tr <- dat[dat$Cohort == tr, , drop = FALSE]
    y_tr   <- factor(dat_tr[[label_col]])
    y_tr   <- factor(y_tr, levels = get_control_first_levels(y_tr, pos_class))
    X_tr   <- dat_tr[, feat_cols, drop = FALSE]

    if (length(unique(y_tr)) < 2) {

      next
    }

    params <- get_representative_params(task_prefix, tr)

    set.seed(100 + which(cohorts == tr))
    rf_tr <- randomForest(
      x        = X_tr,
      y        = y_tr,
      mtry     = params$mtry,
      ntree    = params$ntree,
      nodesize = params$nodesize,
      maxnodes = params$maxnodes,
      importance = TRUE
    )

    for (te in cohorts) {
      if (te == tr) next

      dat_te <- dat[dat$Cohort == te, , drop = FALSE]
      y_te   <- factor(dat_te[[label_col]])
      y_te   <- factor(y_te, levels = levels(y_tr))
      X_te   <- dat_te[, feat_cols, drop = FALSE]

      pair_name <- paste(tr, te, sep = "_")

      if (length(unique(y_te)) < 2) {

        metrics_mat[pair_name, ] <- NA
        next
      }

      pred_class <- predict(rf_tr, X_te)
      pred_prob  <- predict(rf_tr, X_te, type = "prob")[, pos_class]

      cm <- confusionMatrix(
        data      = factor(pred_class, levels = levels(y_tr)),
        reference = y_te,
        mode      = "everything",
        positive  = pos_class
      )

      m_vec <- numeric(length(metric_names))
      m_vec[1:11] <- as.numeric(cm$byClass)

      auc_val <- get_auc(pred_prob, y_te, pos_class)
      m_vec[12] <- auc_val

      metrics_mat[pair_name, ] <- m_vec
    }
  }

  out_file <- file.path(out_dir, paste0(task_prefix, "_CrossCohort_RF_metrics.csv"))
  write.csv(metrics_mat, file = out_file)
}

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

tune_LOCO_model <- function(X, y_factor, pos_class, seed_base = 500) {
  p <- ncol(X)

  set.seed(seed_base)

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
    .ntree    = c(301, 501, 801),
    .nodesize = c(50, 100),
    .maxnodes = c(5, 10, 15)
  )

  y_factor <- factor(y_factor, levels = get_control_first_levels(y_factor, pos_class))

  rf_fit <- caret::train(
    x         = X,
    y         = y_factor,
    method    = customRF,
    metric    = "ROC",
    preProcess = c("center", "scale"),
    tuneGrid  = tunegrid,
    trControl = ctrl
  )

  rf_fit$bestTune
}

run_LOCO <- function(dat, task_prefix, label_col, pos_class, out_dir) {

  feat_cols <- setdiff(colnames(dat), c("Sample","Cohort", label_col))
  cohorts   <- sort(unique(dat$Cohort))

  loco_param_all <- data.frame()
  loco_metrics   <- matrix(NA, nrow = length(cohorts), ncol = length(metric_names))
  rownames(loco_metrics) <- cohorts
  colnames(loco_metrics) <- metric_names

  for (k in cohorts) {

    dat_te <- dat[dat$Cohort == k, , drop = FALSE]
    dat_tr <- dat[dat$Cohort != k, , drop = FALSE]

    y_tr <- factor(dat_tr[[label_col]])
    y_tr <- factor(y_tr, levels = get_control_first_levels(y_tr, pos_class))
    X_tr <- dat_tr[, feat_cols, drop = FALSE]

    y_te <- factor(dat_te[[label_col]])
    y_te <- factor(y_te, levels = levels(y_tr))
    X_te <- dat_te[, feat_cols, drop = FALSE]

    if (length(unique(y_tr)) < 2 || length(unique(y_te)) < 2) {

      next
    }

    best_LOCO <- tune_LOCO_model(
      X         = X_tr,
      y_factor  = y_tr,
      pos_class = pos_class,
      seed_base = 700 + which(cohorts == k)
    )

    loco_param_all <- rbind(
      loco_param_all,
      cbind(Cohort = k, as.data.frame(best_LOCO))
    )

    set.seed(800 + which(cohorts == k))
    rf_LOCO <- randomForest(
      x        = X_tr,
      y        = y_tr,
      mtry     = best_LOCO$mtry,
      ntree    = best_LOCO$ntree,
      nodesize = best_LOCO$nodesize,
      maxnodes = best_LOCO$maxnodes,
      importance = TRUE
    )

    pred_class <- predict(rf_LOCO, X_te)
    pred_prob  <- predict(rf_LOCO, X_te, type = "prob")[, pos_class]

    cm <- confusionMatrix(
      data      = factor(pred_class, levels = levels(y_tr)),
      reference = y_te,
      mode      = "everything",
      positive  = pos_class
    )

    m_vec <- numeric(length(metric_names))
    m_vec[1:11] <- as.numeric(cm$byClass)
    m_vec[12]   <- get_auc(pred_prob, y_te, pos_class)

    loco_metrics[k, ] <- m_vec
  }

  rownames(loco_param_all) <- NULL
  out_param <- file.path(out_dir, paste0(task_prefix, "_LOCO_RF_bestTune.csv"))
  write.csv(loco_param_all, file = out_param, row.names = FALSE)

  out_metric <- file.path(out_dir, paste0(task_prefix, "_LOCO_RF_metrics.csv"))
  write.csv(loco_metrics, file = out_metric)
}

  dat <- read.table(
    TRAIN_DATA,
    header      = TRUE,
    sep         = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  check_cohorts(dat$Cohort)

  run_cross_cohort(
    dat        = dat,
    task_prefix= task_prefix,
    label_col  = label_col,
    pos_class  = pos_class,
    out_dir    = OUT_DIR_CC
  )

  run_LOCO(
    dat        = dat,
    task_prefix= task_prefix,
    label_col  = label_col,
    pos_class  = pos_class,
    out_dir    = OUT_DIR_LOCO
  )
