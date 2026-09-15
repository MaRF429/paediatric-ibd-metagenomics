#!/usr/bin/env Rscript
# Run one comparison and feature type. Input files must already be prepared.
# Usage: Rscript 02_Permutation_manual.R COMPARISON FEATURE_TYPE INPUT_ROOT OUTPUT_ROOT
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
OUT_DIR <- file.path(RESULT_ROOT, "02_Permutation_manual")
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

B_perm <- 1000
get_representative_params <- function(task_prefix, cohort) {
  file_para <- file.path(
    LOCALCV_DIR,
    paste0(task_prefix, "_LocalCV_RF_", cohort, "_selfcv_repeat_bestTune.csv")
  )
  if (!file.exists(file_para)) {
    stop("Required input is unavailable or invalid.", file_para)
  }

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

  res_all <- data.frame()

  for (s in cohorts) {

    dat_s <- dat[dat$Cohort == s, , drop = FALSE]

    y_s <- factor(dat_s[[label_col]])
    y_s <- factor(y_s, levels = get_control_first_levels(y_s, pos_class))

    if (length(unique(y_s)) < 2) {

      next
    }

    X_s <- dat_s[, feat_cols, drop = FALSE]

    params <- get_representative_params(task_prefix = task_prefix, cohort = s)

    set.seed(2021)
    rf_true <- randomForest(
      x        = X_s,
      y        = y_s,
      mtry     = params$mtry,
      ntree    = params$ntree,
      nodesize = params$nodesize,
      maxnodes = params$maxnodes,
      importance = TRUE
    )

    prob_true <- predict(rf_true, X_s, type = "prob")[, pos_class]
    auc_obs   <- get_auc(prob_true, y_s, pos_class)

    auc_perm <- numeric(B_perm)

    set.seed(2022)

    for (b in 1:B_perm) {
      y_perm <- sample(y_s)

      rf_perm <- randomForest(
        x        = X_s,
        y        = y_perm,
        mtry     = params$mtry,
        ntree    = params$ntree,
        nodesize = params$nodesize,
        maxnodes = params$maxnodes
      )

      prob_perm <- predict(rf_perm, X_s, type = "prob")[, pos_class]
      auc_perm[b] <- get_auc(prob_perm, y_perm, pos_class)
    }

    p_perm <- (sum(auc_perm >= auc_obs) + 1) / (B_perm + 1)

    res_s <- data.frame(
      Task       = task_prefix,
      Cohort     = s,
      N_samples  = nrow(dat_s),
      AUC_obs    = auc_obs,
      AUC_perm_mean = mean(auc_perm),
      AUC_perm_sd   = sd(auc_perm),
      p_perm_AUC = p_perm
    )

    res_all <- rbind(res_all, res_s)

    out_null_file <- file.path(
      OUT_DIR,
      paste0(task_prefix, "_Permutation_RF_AUC_null_", s, ".csv")
    )
    write.csv(
      data.frame(AUC_perm = auc_perm),
      file      = out_null_file,
      row.names = FALSE
    )
  }

  out_sum_file <- file.path(
    OUT_DIR,
    paste0(task_prefix, "_Permutation_RF_summary.csv")
  )
  write.csv(
    res_all,
    file      = out_sum_file,
    row.names = FALSE
  )
