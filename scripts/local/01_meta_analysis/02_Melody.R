#!/usr/bin/env Rscript
# Run one comparison and feature type. Input files must already be prepared.
# Usage: Rscript 02_Melody.R COMPARISON FEATURE_TYPE INPUT_ROOT OUTPUT_ROOT
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
OUT_DIR <- file.path(RESULT_ROOT, "Melody")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
pos_class <- if (comparison == "IBD_HC") "IBD" else "CD"
neg_class <- if (comparison == "IBD_HC") "HC" else "UC"

# Enforce the discovery cohort boundary.
check_cohorts <- function(cohorts) {
  allowed <- c("PRJNA398089", "PRJNA1265906", "PRJNA1045596", "PRJNA389280", "HRA007915")
  cohorts <- sub("^(ENA-|GSA-)", "", as.character(cohorts))
  stopifnot(length(cohorts) > 0, !anyNA(cohorts), all(cohorts %in% allowed))
}
library(miMeta)
data <- read.table(
  file.path(INPUT_DIR, "counts.tsv"),
  header = TRUE, row.names = 1, sep = "\t", check.names = FALSE
)
meta <- read.csv(
  file.path(INPUT_DIR, "metadata.csv"),
  row.names = 1)

check_cohorts(meta$Cohort)
data[!is.finite(as.matrix(data))] <- 0

common_samples <- intersect(colnames(data), rownames(meta))
data <- data[, common_samples, drop = FALSE]
meta <- meta[common_samples, , drop = FALSE]

meta$Group <- factor(meta$Group, levels = c("HC", "CD", "UC"))
meta$IBD   <- factor(ifelse(meta$Group %in% c("CD", "UC"), "IBD", "HC"),
                     levels = c("HC", "IBD"))

abd <- t(as.matrix(data))
mode(abd) <- "numeric"

build_melody_input <- function(abd, meta, group_var, case_label, control_label,
                               adjust_vars = c("Age", "Sex")) {
  g <- meta[[group_var]]

  keep <- g %in% c(control_label, case_label) &
    !is.na(g) & !is.na(meta$Cohort)
  meta_sub <- droplevels(meta[keep, , drop = FALSE])
  abd_sub  <- abd[keep, , drop = FALSE]

  nonzero  <- rowSums(abd_sub) > 0
  meta_sub <- meta_sub[nonzero, , drop = FALSE]
  abd_sub  <- abd_sub[nonzero, , drop = FALSE]

  cohorts <- unique(meta_sub$Cohort)
  rel.abd            <- vector("list", length(cohorts))
  covariate.interest <- vector("list", length(cohorts))
  covariates.adjust  <- vector("list", length(cohorts))
  names(rel.abd)            <- cohorts
  names(covariate.interest) <- cohorts
  names(covariates.adjust)  <- cohorts

  for (co in cohorts) {
    idx   <- meta_sub$Cohort == co
    samps <- rownames(meta_sub)[idx]

    m <- as.matrix(abd_sub[idx, , drop = FALSE])
    mode(m) <- "numeric"
    rel.abd[[co]] <- m

    grp     <- meta_sub[[group_var]][idx]
    disease <- ifelse(grp == case_label, 1, 0)
    disease <- as.numeric(disease)
    names(disease) <- samps
    covariate.interest[[co]] <- matrix(
      disease,
      ncol = 1,
      dimnames = list(samps, "disease")
    )

    adj <- meta_sub[idx, adjust_vars, drop = FALSE]

    adj$Age <- as.numeric(as.character(adj$Age))
    adj$Sex <- ifelse(adj$Sex %in% c("Male", "M", 1), 1,
                      ifelse(adj$Sex %in% c("Female", "F", 0), 0, NA))
    adj$Sex <- as.numeric(adj$Sex)

    rownames(adj) <- samps
    covariates.adjust[[co]] <- adj
  }

  list(
    rel.abd            = rel.abd,
    covariate.interest = covariate.interest,
    covariates.adjust  = covariates.adjust
  )
}

clean_covariates <- function(input) {
  for (co in names(input$covariates.adjust)) {
    adj <- input$covariates.adjust[[co]]
    bad_cols <- vapply(adj, function(x) any(is.na(x)), logical(1))
    adj <- adj[, !bad_cols, drop = FALSE]

    if (ncol(adj) == 0) {
      input$covariates.adjust[[co]] <- NULL
      next
    }

    input$covariates.adjust[[co]] <- adj
  }
  return(input)
}

filter_studies_by_group_n <- function(input, min_n = 5) {
  keep <- names(input$covariate.interest)[
    vapply(input$covariate.interest, function(Z) {
      z <- Z[, "disease"]
      n_case <- sum(z == 1, na.rm = TRUE)
      n_ctrl <- sum(z == 0, na.rm = TRUE)
      (n_case > min_n) && (n_ctrl > min_n)
    }, logical(1))
  ]
  input$rel.abd            <- input$rel.abd[keep]
  input$covariate.interest <- input$covariate.interest[keep]
  input$covariates.adjust  <- input$covariates.adjust[keep]
  return(input)
}

input <- build_melody_input(abd, meta,
  group_var = if (comparison == "IBD_HC") "IBD" else "Group",
  case_label = pos_class, control_label = neg_class, adjust_vars = c("Age", "Sex"))
input <- clean_covariates(input)
input <- filter_studies_by_group_n(input, min_n = 5)
fit <- melody(rel.abd = input$rel.abd,
  covariate.interest = input$covariate.interest,
  covariate.adjust = input$covariates.adjust,
  prev.filter = 0.1, ref = NULL, verbose = TRUE)
null <- melody.null.model(rel.abd = input$rel.abd,
  covariate.adjust = input$covariates.adjust, prev.filter = 0.1, ref = NULL)
summary <- melody.get.summary(null.obj = null,
  covariate.interest = input$covariate.interest, cluster = NULL, verbose = TRUE)
coefficients <- fit$disease$coef
result <- data.frame(Feature = names(coefficients),
  Melody_coef = as.numeric(coefficients),
  Melody_selected = as.integer(coefficients != 0))
write.csv(result, file.path(OUT_DIR, "Melody_coefficients.csv"), row.names = FALSE)
saveRDS(summary, file.path(OUT_DIR, "Melody_summary.rds"))
