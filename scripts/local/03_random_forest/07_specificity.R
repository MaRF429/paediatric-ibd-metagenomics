#!/usr/bin/env Rscript
# Run one comparison and feature type. Input files must already be prepared.
# Usage: Rscript 07_specificity.R IBD_HC FEATURE_TYPE INPUT_ROOT OUTPUT_ROOT SHARED_CONTROL_COHORT
# Inputs and outputs are isolated under COMPARISON/FEATURE_TYPE.
options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("Expected comparison, feature type, input root and output root.")
comparison <- args[[1]]
feature_type <- args[[2]]
stopifnot(comparison == "IBD_HC",
          feature_type %in% c("species", "KO", "EggNOG", "pathway", "combined"))

INPUT_DIR <- file.path(args[[3]], comparison, feature_type)
RESULT_ROOT <- file.path(args[[4]], comparison, feature_type)
OUT_DIR <- file.path(RESULT_ROOT, "07_specificity")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
task_prefix <- "IBD"
TEST_DATA <- file.path(INPUT_DIR, paste0(task_prefix, "_test_data_use.txt"))
MODEL_PATH <- file.path(RESULT_ROOT, "04_FinalModel", paste0(task_prefix, "_Final_RF_model.rds"))

# Enforce the independent-validation cohort boundary.
check_cohorts <- function(cohorts) {
  allowed <- c("PRJEB76677", "SRP057027", "PRJNA759642", "PRJNA922068")
  cohorts <- sub("^(ENA-|GSA-)", "", as.character(cohorts))
  stopifnot(length(cohorts) > 0, !anyNA(cohorts), all(cohorts %in% allowed))
}

# Fifth argument names the study's shared-control challenge cohort.
if (length(args) < 5 || !nzchar(args[[5]])) stop("Specify the shared-control cohort as argument 5.")
library(caret)
library(randomForest)
library(ROCR)
set.seed(123)
SELF_TEST_COHORT <- args[[5]]
MODEL_POS_CLASS <- "IBD"
NEG_CLASS <- "HC"
GROUPS_TO_EVALUATE <- c("ASD", "DM", "KD", "SID", "IBD")

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

get_control_first_levels <- function(y_factor, pos_class) {
  c(setdiff(levels(factor(y_factor)), pos_class)[1], pos_class)
}

evaluate_binary <- function(dat, pos_class) {
  y <- factor(dat$Eval_Label)
  y <- factor(y, levels = get_control_first_levels(y, pos_class))
  pred_class <- factor(dat$Eval_Pred, levels = levels(y))

  cm <- confusionMatrix(
    data = pred_class,
    reference = y,
    positive = pos_class,
    mode = "everything"
  )

  data.frame(
    Sensitivity = unname(cm$byClass["Sensitivity"]),
    Specificity = unname(cm$byClass["Specificity"]),
    `Pos Pred Value` = unname(cm$byClass["Pos Pred Value"]),
    `Neg Pred Value` = unname(cm$byClass["Neg Pred Value"]),
    Precision = unname(cm$byClass["Precision"]),
    Recall = unname(cm$byClass["Recall"]),
    F1 = unname(cm$byClass["F1"]),
    Prevalence = unname(cm$byClass["Prevalence"]),
    `Detection Rate` = unname(cm$byClass["Detection Rate"]),
    `Detection Prevalence` = unname(cm$byClass["Detection Prevalence"]),
    `Balanced Accuracy` = unname(cm$byClass["Balanced Accuracy"]),
    AUC = get_auc(dat$Prob_IBD, y, pos_class),
    check.names = FALSE
  )
}

get_reference_hc <- function(pred_df, group_name) {
  if (group_name %in% c("DM", "KD", "SID")) {
    return(pred_df[pred_df$Cohort == SELF_TEST_COHORT & pred_df$Group == NEG_CLASS, , drop = FALSE])
  }

  cohort_use <- unique(pred_df$Cohort[pred_df$Group == group_name])
  pred_df[pred_df$Cohort %in% cohort_use & pred_df$Group == NEG_CLASS, , drop = FALSE]
}

evaluate_group <- function(pred_df, group_name) {
  group_dat <- pred_df[pred_df$Group == group_name, , drop = FALSE]
  ref_hc <- get_reference_hc(pred_df, group_name)

  if (nrow(group_dat) == 0) stop("No samples found for group: ", group_name)
  if (nrow(ref_hc) == 0) stop("No matched HC found for group: ", group_name)

  eval_dat <- rbind(group_dat, ref_hc)

  eval_dat$Eval_Label <- ifelse(eval_dat$Group == group_name, group_name, NEG_CLASS)
  eval_dat$Eval_Pred <- ifelse(eval_dat$Pred_Class == MODEL_POS_CLASS, group_name, NEG_CLASS)

  metrics <- evaluate_binary(eval_dat, group_name)

  data.frame(
    Group = group_name,
    Validation_Cohort = paste(sort(unique(eval_dat$Cohort)), collapse = ";"),
    Reference_HC = paste(sort(unique(ref_hc$Cohort)), collapse = ";"),
    Positive_Label = group_name,
    Group_N = nrow(group_dat),
    Reference_HC_N = nrow(ref_hc),
    Group_Pred_IBD_N = sum(group_dat$Pred_Class == MODEL_POS_CLASS),
    Group_Pred_IBD_rate = mean(group_dat$Pred_Class == MODEL_POS_CLASS),
    Group_Mean_Prob_IBD = mean(group_dat$Prob_IBD),
    Group_Median_Prob_IBD = median(group_dat$Prob_IBD),
    Reference_HC_Pred_IBD_rate = mean(ref_hc$Pred_Class == MODEL_POS_CLASS),
    metrics,
    check.names = FALSE
  )
}

dat <- read.table(file.path(INPUT_DIR, "specificity_data.tsv"), header = TRUE, sep = "\t", check.names = FALSE)
rf_model <- readRDS(MODEL_PATH)

feat_cols <- names(rf_model$forest$xlevels)
X <- dat[, feat_cols, drop = FALSE]

prob_mat <- predict(rf_model, X, type = "prob")
prob_pos <- prob_mat[, MODEL_POS_CLASS]

pred_df <- data.frame(
  Sample = dat$Sample,
  Cohort = dat$Cohort,
  Group = dat$Group,
  Group_IBD = dat$Group_IBD,
  Design = ifelse(dat$Cohort == SELF_TEST_COHORT,
                  "Self_test_with_shared_HC",
                  "External_cohort_with_own_HC"),
  Pred_Class = factor(ifelse(prob_pos >= 0.5, MODEL_POS_CLASS, NEG_CLASS),
                      levels = c(NEG_CLASS, MODEL_POS_CLASS)),
  Prob_HC = prob_mat[, NEG_CLASS],
  Prob_IBD = prob_pos,
  stringsAsFactors = FALSE
)

metrics_df <- do.call(rbind, lapply(GROUPS_TO_EVALUATE, function(g) evaluate_group(pred_df, g)))
rownames(metrics_df) <- NULL

write.csv(
  pred_df,
  file = file.path(OUT_DIR, "IBD_specificity_FinalRF_predictions.csv"),
  row.names = FALSE
)

write.csv(
  metrics_df,
  file = file.path(OUT_DIR, "IBD_specificity_FinalRF_metrics.csv"),
  row.names = FALSE
)

# Expanded-control AUROC; reset the original seed for this second stage.
set.seed(123)
ADD_N <- 5
N_REPEAT <- 50

challenge_sets <- data.frame(
  Disease = c("ASD", "DM", "KD", "SID", "adult-IBD"),
  Challenge = c("ASD challenge", "DM challenge", "KD challenge",
                "SID challenge", "adult-IBD challenge"),
  Patient_Cohort = c("OEP00004172", SELF_TEST_COHORT, SELF_TEST_COHORT,
                     SELF_TEST_COHORT, "PRJCA017408"),
  Patient_Group = c("ASD", "DM", "KD", "SID", "IBD"),
  HC_Cohort = c("OEP00004172", SELF_TEST_COHORT, SELF_TEST_COHORT,
                SELF_TEST_COHORT, "PRJCA017408"),
  HC_Group = rep("HC", 5),
  stringsAsFactors = FALSE
)

disease_levels <- c("DM", "KD", "SID", "ASD", "adult-IBD")

predict_prob_ibd <- function(rf_model, dat, feat_cols) {
  predict(rf_model, dat[, feat_cols, drop = FALSE], type = "prob")[, MODEL_POS_CLASS]
}

make_auc_data <- function(base_pos, base_ctrl, added_ctrl) {

  dat_auc <- rbind(
    data.frame(Prob_IBD = base_pos$Prob_IBD, Label = MODEL_POS_CLASS),
    data.frame(Prob_IBD = base_ctrl$Prob_IBD, Label = NEG_CLASS),
    data.frame(Prob_IBD = added_ctrl$Prob_IBD, Label = rep(NEG_CLASS, nrow(added_ctrl)))
  )
  dat_auc$Label <- factor(dat_auc$Label, levels = c(NEG_CLASS, MODEL_POS_CLASS))
  dat_auc
}

calc_expanded_auc <- function(base_pos, base_ctrl, added_ctrl) {
  auc_dat <- make_auc_data(base_pos, base_ctrl, added_ctrl)
  get_auc(auc_dat$Prob_IBD, auc_dat$Label, MODEL_POS_CLASS)
}

ped_val <- read.table(TEST_DATA, header = TRUE, sep = "\t", check.names = FALSE)
spec_dat <- dat

check_cohorts(ped_val$Cohort)
ped_val$Prob_IBD <- predict_prob_ibd(rf_model, ped_val, feat_cols)
spec_dat$Prob_IBD <- predict_prob_ibd(rf_model, spec_dat, feat_cols)

base_pos <- ped_val[ped_val$Group_IBD == MODEL_POS_CLASS, , drop = FALSE]
base_ctrl <- ped_val[ped_val$Group_IBD == NEG_CLASS, , drop = FALSE]

original_auc <- calc_expanded_auc(base_pos, base_ctrl, base_ctrl[0, , drop = FALSE])

paired_list <- list()

for (i in seq_len(nrow(challenge_sets))) {
  disease <- challenge_sets$Disease[i]
  challenge <- challenge_sets$Challenge[i]

  patient_pool <- spec_dat[
    spec_dat$Cohort == challenge_sets$Patient_Cohort[i] &
      spec_dat$Group == challenge_sets$Patient_Group[i],
    ,
    drop = FALSE
  ]

  hc_pool <- spec_dat[
    spec_dat$Cohort == challenge_sets$HC_Cohort[i] &
      spec_dat$Group == challenge_sets$HC_Group[i],
    ,
    drop = FALSE
  ]

  if (nrow(patient_pool) < ADD_N) {
    stop(disease, " patient pool has fewer than ADD_N samples.")
  }
  if (nrow(hc_pool) < ADD_N) {
    stop(disease, " HC pool has fewer than ADD_N samples.")
  }

  patient_auc <- numeric(N_REPEAT)
  control_auc <- numeric(N_REPEAT)

  for (r in seq_len(N_REPEAT)) {
    patient_idx <- sample(seq_len(nrow(patient_pool)), size = ADD_N, replace = FALSE)
    hc_idx <- sample(seq_len(nrow(hc_pool)), size = ADD_N, replace = FALSE)

    patient_auc[r] <- calc_expanded_auc(base_pos, base_ctrl, patient_pool[patient_idx, , drop = FALSE])
    control_auc[r] <- calc_expanded_auc(base_pos, base_ctrl, hc_pool[hc_idx, , drop = FALSE])
  }

  paired_list[[disease]] <- data.frame(
    Disease = disease,
    Challenge = challenge,
    Patient_Cohort = challenge_sets$Patient_Cohort[i],
    Patient_Group = challenge_sets$Patient_Group[i],
    HC_Cohort = challenge_sets$HC_Cohort[i],
    HC_Group = challenge_sets$HC_Group[i],
    Patient_Pool_N = nrow(patient_pool),
    HC_Pool_N = nrow(hc_pool),
    Add_N = ADD_N,
    Repeat = seq_len(N_REPEAT),
    Original_AUC = original_auc,
    Patient_AUC = patient_auc,
    Control_AUC = control_auc,
    Delta_AUC = patient_auc - control_auc,
    Delta_percent = 100 * (patient_auc - control_auc) / control_auc,
    Added_Label = NEG_CLASS,
    Model_Retrained = "No",
    stringsAsFactors = FALSE
  )
}

paired_df <- do.call(rbind, paired_list)
rownames(paired_df) <- NULL
paired_df$Disease <- factor(paired_df$Disease, levels = disease_levels)

summary_df <- do.call(rbind, lapply(split(paired_df, paired_df$Disease), function(x) {
  data.frame(
    Disease = as.character(unique(x$Disease)),
    Challenge = unique(x$Challenge),
    Patient_Cohort = unique(x$Patient_Cohort),
    Patient_Group = unique(x$Patient_Group),
    HC_Cohort = unique(x$HC_Cohort),
    HC_Group = unique(x$HC_Group),
    Patient_Pool_N = unique(x$Patient_Pool_N),
    HC_Pool_N = unique(x$HC_Pool_N),
    Add_N = unique(x$Add_N),
    N_REPEAT = nrow(x),
    Original_AUC = unique(x$Original_AUC),
    Mean_Patient_AUC = mean(x$Patient_AUC),
    Median_Patient_AUC = median(x$Patient_AUC),
    Mean_Control_AUC = mean(x$Control_AUC),
    Median_Control_AUC = median(x$Control_AUC),
    Mean_Delta_AUC = mean(x$Delta_AUC),
    Median_Delta_AUC = median(x$Delta_AUC),
    Mean_Delta_percent = mean(x$Delta_percent),
    Median_Delta_percent = median(x$Delta_percent),
    Added_Label = unique(x$Added_Label),
    Model_Retrained = "No",
    stringsAsFactors = FALSE
  )
}))
rownames(summary_df) <- NULL
summary_df$Disease <- factor(summary_df$Disease, levels = disease_levels)
summary_df <- summary_df[order(summary_df$Disease), ]

write.csv(
  paired_df,
  file = file.path(OUT_DIR, "IBD_expanded_control_AUC_repeats.csv"),
  row.names = FALSE
)

write.csv(
  summary_df,
  file = file.path(OUT_DIR, "IBD_expanded_control_AUC_summary.csv"),
  row.names = FALSE
)
