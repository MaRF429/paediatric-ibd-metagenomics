#!/usr/bin/env Rscript
# Run one comparison and feature type. Input files must already be prepared.
# Usage: Rscript 01_MMUPHin.R COMPARISON FEATURE_TYPE INPUT_ROOT OUTPUT_ROOT
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
OUT_DIR <- file.path(RESULT_ROOT, "MMUPHin")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
pos_class <- if (comparison == "IBD_HC") "IBD" else "CD"
neg_class <- if (comparison == "IBD_HC") "HC" else "UC"

# Enforce the discovery cohort boundary.
check_cohorts <- function(cohorts) {
  allowed <- c("PRJNA398089", "PRJNA1265906", "PRJNA1045596", "PRJNA389280", "HRA007915")
  cohorts <- sub("^(ENA-|GSA-)", "", as.character(cohorts))
  stopifnot(length(cohorts) > 0, !anyNA(cohorts), all(cohorts %in% allowed))
}
library(MMUPHin)
library(dplyr)
library(tidyr)
data <- read.table(
  file.path(INPUT_DIR, "abundance.tsv"),
  header = TRUE, row.names = 1, sep = "\t", check.names = FALSE
)
meta <- read.csv(
  file.path(INPUT_DIR, "metadata.csv"),
  row.names = 1)

data <- as.matrix(data)
stopifnot(is.numeric(data[1, 1]) || is.integer(data[1, 1]))
stopifnot(all(data >= 0 & data <= 1, na.rm = TRUE))

check_cohorts(meta$Cohort)
meta$Cohort <- factor(meta$Cohort)
meta$Sex    <- factor(ifelse(is.na(meta$Sex) | meta$Sex == "", "Unknown", meta$Sex))
meta$Age    <- suppressWarnings(as.numeric(meta$Age))
meta$Group  <- factor(meta$Group, levels = c("HC", "CD", "UC"))
meta$IBD    <- factor(ifelse(meta$Group %in% c("CD", "UC"), "IBD", "HC"),
                      levels = c("HC", "IBD"))

common_ids <- intersect(colnames(data), rownames(meta))
stopifnot(length(common_ids) > 0)
data  <- data[,  common_ids, drop = FALSE]
meta2 <- meta[common_ids, , drop = FALSE]

meta2$IBD   <- stats::relevel(meta2$IBD,   ref = "HC")
meta2$Group <- stats::relevel(meta2$Group, ref = "HC")

run_contrast <- function(feature_abd, sample_anno,
                         exposure, keep_levels, ref,
                         tag, out_dir,
                         min_per_group = 5) {

  stopifnot(exposure %in% colnames(sample_anno),
            "Cohort" %in% colnames(sample_anno))

  ids <- rownames(sample_anno)[sample_anno[[exposure]] %in% keep_levels]
  ids <- intersect(ids, colnames(feature_abd))
  stopifnot(length(ids) > 0)

  abd2  <- feature_abd[, ids, drop = FALSE]
  anno2 <- sample_anno[ids, , drop = FALSE]

  anno2[[exposure]] <- as.character(anno2[[exposure]])
  anno2$Cohort       <- as.character(anno2$Cohort)

  ok <- !is.na(anno2[[exposure]]) & !is.na(anno2$Cohort) & (anno2[[exposure]] %in% keep_levels)
  anno2 <- anno2[ok, , drop = FALSE]
  abd2  <- abd2[, rownames(anno2), drop = FALSE]
  stopifnot(nrow(anno2) > 0)

  tab <- anno2 |>
    dplyr::count(Cohort, grp = .data[[exposure]]) |>
    tidyr::pivot_wider(names_from = grp, values_from = n, values_fill = 0)

  for (lv in keep_levels) if (!lv %in% names(tab)) tab[[lv]] <- 0L

  keep_cohorts <- tab |>
    dplyr::filter(.data[[keep_levels[1]]] > min_per_group,
                  .data[[keep_levels[2]]] > min_per_group) |>
    dplyr::pull(Cohort)

  anno2 <- anno2[anno2$Cohort %in% keep_cohorts, , drop = FALSE]
  abd2  <- abd2[, rownames(anno2), drop = FALSE]
  if (nrow(anno2) == 0) {
    stop("No cohorts left after cohort filter (both groups > ",
         min_per_group, ") for tag=", tag)
  }

  anno2[[exposure]] <- factor(anno2[[exposure]], levels = keep_levels)

    anno2[[exposure]] <- stats::relevel(anno2[[exposure]], ref = ref)

  anno2$Age <- suppressWarnings(as.numeric(anno2$Age))

    anno2$Sex <- as.character(anno2$Sex)
    anno2$Sex[is.na(anno2$Sex) | anno2$Sex == ""] <- "Unknown"
    anno2$Sex <- factor(anno2$Sex)

  tab2 <- table(anno2[[exposure]])
  stopifnot(all(tab2 > 0))

  fit <- MMUPHin::lm_meta(
    feature_abd = abd2,
    batch       = "Cohort",
    exposure    = exposure,
    covariates  = c("Sex", "Age"),
    data        = anno2,
    control     = list(
      normalization   = "NONE",
      transform       = "AST",
      analysis_method = "LM",
      rma_method      = "REML",
      verbose         = TRUE
    )
  )

  res <- fit$meta_fits[order(fit$meta_fits$qval.fdr, fit$meta_fits$pval), ]
  write.csv(res, file.path(out_dir, paste0("mmuphin_", tag, ".csv")), row.names = FALSE)

  invisible(list(res = res, fit = fit,
                 kept_cohorts = unique(anno2$Cohort),
                 n_samples = ncol(abd2)))
}

gather_maaslin_meta <- function(fit, exposure) {
  cohort_long <- dplyr::bind_rows(fit$maaslin_fits, .id = "cohort") %>%
    dplyr::filter(metadata == exposure) %>%
    dplyr::select(feature, cohort, coef, pval, qval)

  by_cohort <- tidyr::pivot_wider(
    cohort_long,
    id_cols    = feature,
    names_from = cohort,
    values_from = c(coef, pval, qval),
    names_sep = "_"
  )

  meta_out <- fit$meta_fits %>%
    dplyr::select(
      feature,
      meta_coef = coef,
      meta_pval = pval,
      meta_qval = qval.fdr,
      meta_I2   = I2
    )

  meta_out %>%
    dplyr::left_join(by_cohort, by = "feature") %>%
    dplyr::arrange(meta_qval, meta_pval)
}

filter_meta_direction <- function(tbl,
                                  coef_prefix = "coef_",
                                  alpha = 0.05,
                                  min_agree = 4,
                                  require_meta_match = FALSE) {
  coef_cols <- grep(paste0("^", coef_prefix), names(tbl), value = TRUE)
  stopifnot(length(coef_cols) > 0,
            all(c("meta_coef", "meta_qval", "meta_pval") %in% names(tbl)))

  tbl %>%
    dplyr::mutate(
      pos = rowSums(dplyr::across(all_of(coef_cols), ~ .x > 0), na.rm = TRUE),
      neg = rowSums(dplyr::across(all_of(coef_cols), ~ .x < 0), na.rm = TRUE),
      agree_n       = pmax(pos, neg),
      majority_sign = sign(pos - neg),
      cond_q        = meta_qval < alpha,
      cond_dir      = agree_n   >= min_agree,
      cond_meta     = (!require_meta_match) | (majority_sign == sign(meta_coef)),
      pass          = cond_q & cond_dir & cond_meta
    ) %>%
    dplyr::filter(pass) %>%
    dplyr::arrange(meta_qval, meta_pval)
}

exposure <- if (comparison == "IBD_HC") "IBD" else "Group"
min_agree <- if (comparison == "IBD_HC") 5 else 3
fit <- run_contrast(data, meta2, exposure, c(neg_class, pos_class),
                    neg_class, comparison, OUT_DIR, min_per_group = 5)
tbl <- gather_maaslin_meta(fit$fit, exposure)
all_results <- filter_meta_direction(tbl, alpha = 1, min_agree = 0,
                                     require_meta_match = TRUE)
selected <- filter_meta_direction(tbl, alpha = 0.05, min_agree = min_agree,
                                  require_meta_match = TRUE)
readr::write_csv(all_results, file.path(OUT_DIR, "all_results.csv"))
readr::write_csv(selected, file.path(OUT_DIR, "filtered_results.csv"))
