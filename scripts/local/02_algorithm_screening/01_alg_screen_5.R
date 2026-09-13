#!/usr/bin/env Rscript
# Run one comparison and feature type. Input files must already be prepared.
# Usage: Rscript 01_alg_screen_5.R COMPARISON FEATURE_TYPE INPUT_ROOT OUTPUT_ROOT
# Inputs and outputs are isolated under COMPARISON/FEATURE_TYPE.
options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("Expected comparison, feature type, input root and output root.")
comparison <- args[[1]]
feature_type <- args[[2]]
stopifnot(comparison %in% c("IBD_HC", "CD_UC"),
          feature_type %in% c("species", "KO", "EggNOG", "pathway"))

INPUT_DIR <- file.path(args[[3]], comparison, feature_type)
RESULT_ROOT <- file.path(args[[4]], comparison, feature_type)
OUT_DIR <- file.path(RESULT_ROOT, "algorithm_screening_5")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
pos_class <- if (comparison == "IBD_HC") "IBD" else "CD"
neg_class <- if (comparison == "IBD_HC") "HC" else "UC"

# Enforce the discovery cohort boundary.
check_cohorts <- function(cohorts) {
  allowed <- c("PRJNA398089", "PRJNA1265906", "PRJNA1045596", "PRJNA389280", "HRA007915")
  cohorts <- sub("^(ENA-|GSA-)", "", as.character(cohorts))
  stopifnot(length(cohorts) > 0, !anyNA(cohorts), all(cohorts %in% allowed))
}

library(caret)
library(dplyr)
library(tidyr)
set.seed(123)
# Panel preprocessing and factor levels must match the original comparison.
panel <- readRDS(file.path(INPUT_DIR, "panel.rds"))
metadata <- read.csv(file.path(INPUT_DIR, "metadata.csv"), row.names = 1)
check_cohorts(metadata$Cohort)
X <- panel$X_model
y <- panel$y
stopifnot(all(rownames(X) %in% rownames(metadata)), nrow(X) == length(y),
          all(as.character(y) %in% c(neg_class, pos_class)))
# Preserve the input factor order used by caret::twoClassSummary.
model_settings <- data.frame(
  AlgorithmName  = c("RandomForest", "GradientBoosting", "Lasso",
                     "PLSModel", "DecisionTree"),
  Implementation = c("rf", "xgbTree", "glmnet",
                     "pls", "rpart"),
  stringsAsFactors = FALSE
)

ctrl <- trainControl(
  method          = "repeatedcv",
  number          = 5,
  repeats         = 3,
  savePredictions = "final",
  classProbs      = TRUE,
  summaryFunction = twoClassSummary,
  verboseIter     = FALSE
)

run_models <- function(X, y, model_settings, ctrl,
                       top_features_per_model = 20) {

  model_fits       <- list()
  gene_lists       <- list()
  feature_imp_long <- tibble::tibble()

  for (i in seq_len(nrow(model_settings))) {
    method     <- model_settings$Implementation[i]
    model_name <- model_settings$AlgorithmName[i]

    fit <- caret::train(
      x         = X,
      y         = as.factor(y),
      method    = method,
      trControl = ctrl,
      metric    = "ROC"
    )

    model_fits[[model_name]] <- fit

    imp_raw <- varImp(fit)$importance
    imp <- tibble::tibble(
      Feature    = rownames(imp_raw),
      Importance = as.numeric(imp_raw[, 1]),
      Model      = model_name
    )

    feature_imp_long <- dplyr::bind_rows(feature_imp_long, imp)

    top_feats <- imp %>%
      arrange(desc(Importance)) %>%
      slice_head(n = top_features_per_model) %>%
      pull(Feature)

    gene_lists[[model_name]] <- top_feats
  }

  list(

    model_fits       = model_fits,
    gene_lists       = gene_lists,
    feature_imp_long = feature_imp_long
  )
}

top_features_per_model <- 20
min_model_count_panel <- 3
final_panel_n <- 20
res <- run_models(X, y, model_settings, ctrl, top_features_per_model)
write.csv(res$feature_imp_long, file.path(OUT_DIR, "feature_importance_by_algorithm.csv"),
          row.names = FALSE)
res$feature_imp_long <- res$feature_imp_long %>%
  dplyr::filter(!grepl("^`.*`$", Feature))

wide_imp <- res$feature_imp_long %>%
  dplyr::select(Feature, Model, Importance) %>%
  tidyr::pivot_wider(
    names_from  = Model,
    values_from = Importance
  )

model_cols <- setdiff(names(wide_imp), "Feature")

model_count <- wide_imp %>%
  mutate(
    ModelCount = rowSums(
      dplyr::across(
        all_of(model_cols),
        ~ !is.na(.) & . > 0
      ),
      na.rm = TRUE
    )
  ) %>%
  dplyr::select(Feature, ModelCount)

wide_imp_norm <- wide_imp %>%
  mutate(
    dplyr::across(
      all_of(model_cols),
      ~ {
        rng <- range(., na.rm = TRUE)
        if (!is.finite(rng[1]) || !is.finite(rng[2]) || diff(rng) == 0) {
          rep(0, length(.))
        } else {
          (.- rng[1]) / diff(rng) * 100
        }
      }
    )
  )

feature_summary <- wide_imp_norm %>%
  dplyr::left_join(model_count, by = "Feature") %>%
  mutate(
    MeanImportance = rowMeans(
      dplyr::across(all_of(model_cols)),
      na.rm = TRUE
    )
  ) %>%
  arrange(desc(MeanImportance))

write.csv(feature_summary,
          file.path(OUT_DIR, "feature_importance_summary.csv"),
          row.names = FALSE)

feature_summary_filt <- feature_summary %>%
  filter(ModelCount >= min_model_count_panel)

final_panel <- feature_summary_filt %>%
  slice_head(n = final_panel_n)

write.csv(final_panel,
          file.path(OUT_DIR, "final_panel.csv"),
          row.names = FALSE)

saveRDS(
  list(
    Features        = final_panel$Feature,
    PanelSummary    = final_panel,
    FeatureSummary  = feature_summary,
    ModelSettings   = model_settings
  ),
  file = file.path(OUT_DIR, "final_panel.rds")
)
