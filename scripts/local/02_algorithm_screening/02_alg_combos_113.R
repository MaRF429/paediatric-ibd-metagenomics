#!/usr/bin/env Rscript
# Run one comparison and feature type. Input files must already be prepared.
# Usage: Rscript 02_alg_combos_113.R COMPARISON FEATURE_TYPE INPUT_ROOT OUTPUT_ROOT
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
OUT_DIR <- file.path(RESULT_ROOT, "algorithm_combinations_113")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Enforce the discovery cohort boundary.
check_cohorts <- function(cohorts) {
  allowed <- c("PRJNA398089", "PRJNA1265906", "PRJNA1045596", "PRJNA389280", "HRA007915")
  cohorts <- sub("^(ENA-|GSA-)", "", as.character(cohorts))
  stopifnot(length(cohorts) > 0, !anyNA(cohorts), all(cohorts %in% allowed))
}

library(openxlsx)
library(randomForestSRC)
library(glmnet)
library(plsRglm)
library(gbm)
library(caret)
library(mboost)
library(e1071)
library(MASS)
library(xgboost)
library(pROC)
# Numeric outcome coding is preserved from the supplied comparison-specific input.
Train_expr <- read.table(file.path(INPUT_DIR, "Training_expr.txt"),
  header = TRUE, sep = "\t", row.names = 1, check.names = FALSE,
  stringsAsFactors = FALSE)
Train_class <- read.table(file.path(INPUT_DIR, "Training_class.txt"),
  header = TRUE, sep = "\t", row.names = 1, check.names = FALSE,
  stringsAsFactors = FALSE)
check_cohorts(Train_class$Cohort)
stopifnot(identical(rownames(Train_expr), rownames(Train_class)),
          all(Train_class$outcome %in% c(0, 1)))
# Feature naming is determined exclusively from discovery inputs.
all_feat <- colnames(Train_expr)
safe_feat <- make.names(all_feat, unique = TRUE)
feat_map <- setNames(safe_feat, all_feat)
colnames(Train_expr) <- feat_map[colnames(Train_expr)]
Train_set <- as.matrix(Train_expr)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- dirname(sub("^--file=", "", script_arg[[1]]))
methods <- read.xlsx(file.path(script_dir, "methods.xlsx"), startRow = 2)$Model
methods <- gsub("-| ", "", methods)
stopifnot(length(methods) == 113L, !anyNA(methods), !anyDuplicated(methods))
RunML <- function(method, Train_set, Train_label, mode = "Model", classVar){

  method = gsub(" ", "", method)
  method_name = gsub("(\\w+)\\[(.+)\\]", "\\1", method)
  method_param = gsub("(\\w+)\\[(.+)\\]", "\\2", method)

  method_param = switch(
    EXPR = method_name,
    "Enet" = list("alpha" = as.numeric(gsub("alpha=", "", method_param))),
    "Stepglm" = list("direction" = method_param),
    NULL
  )
  message("Run ", method_name, " algorithm for ", mode, "; ",
          method_param, ";",
          " using ", ncol(Train_set), " Variables")

  args = list("Train_set" = Train_set,
              "Train_label" = Train_label,
              "classVar" = classVar)
  args = c(args, method_param)

  obj <- do.call(what = paste0("Run", method_name),
                 args = args)

  if (mode == "Variable") return(ExtractVar(obj))
  return(obj)
}

RunEnet <- function(Train_set, Train_label, classVar, alpha){
  cv.fit = cv.glmnet(x = Train_set,
                     y = Train_label[[classVar]],
                     family = "binomial", alpha = alpha, nfolds = 10)
  fit = glmnet(x = Train_set,
               y = Train_label[[classVar]],
               family = "binomial", alpha = alpha, lambda = cv.fit$lambda.min)
  fit$subFeature = colnames(Train_set)
  return(fit)
}

RunLasso <- function(Train_set, Train_label, classVar){
  RunEnet(Train_set, Train_label, classVar, alpha = 1)
}

RunRidge <- function(Train_set, Train_label, classVar){
  RunEnet(Train_set, Train_label, classVar, alpha = 0)
}

RunStepglm <- function(Train_set, Train_label, classVar, direction){
  fit <- step(glm(formula = Train_label[[classVar]] ~ .,
                  family = "binomial",
                  data = as.data.frame(Train_set)),
              direction = direction, trace = 0)
  fit$subFeature = colnames(Train_set)
  return(fit)
}

RunSVM <- function(Train_set, Train_label, classVar){
  data <- as.data.frame(Train_set)
  data[[classVar]] <- as.factor(Train_label[[classVar]])
  fit = svm(formula = eval(parse(text = paste(classVar, "~."))),
            data= data, probability = T)
  fit$subFeature = colnames(Train_set)
  return(fit)
}

RunLDA <- function(Train_set, Train_label, classVar){
  data <- as.data.frame(Train_set)
  data[[classVar]] <- as.factor(Train_label[[classVar]])
  fit = train(eval(parse(text = paste(classVar, "~."))),
              data = data,
              method="lda",
              trControl = trainControl(method = "cv"))
  fit$subFeature = colnames(Train_set)
  return(fit)
}

RunglmBoost <- function(Train_set, Train_label, classVar){
  data <- cbind(Train_set, Train_label[classVar])
  data[[classVar]] <- as.factor(data[[classVar]])
  fit <- glmboost(eval(parse(text = paste(classVar, "~."))),
                  data = data,
                  family = Binomial())

  cvm <- cvrisk(fit, papply = lapply,
                folds = cv(model.weights(fit), type = "kfold"))
  fit <- glmboost(eval(parse(text = paste(classVar, "~."))),
                  data = data,
                  family = Binomial(),
                  control = boost_control(mstop = max(mstop(cvm), 40)))

  fit$subFeature = colnames(Train_set)
  return(fit)
}

RunplsRglm <- function(Train_set, Train_label, classVar){
  cv.plsRglm.res = cv.plsRglm(formula = Train_label[[classVar]] ~ .,
                              data = as.data.frame(Train_set),
                              nt=10, verbose = FALSE)
  fit <- plsRglm(Train_label[[classVar]],
                 as.data.frame(Train_set),
                 modele = "pls-glm-logistic",
                 verbose = F, sparse = T)
  fit$subFeature = colnames(Train_set)
  return(fit)
}

RunRF <- function(Train_set, Train_label, classVar){
  rf_nodesize = 5
  Train_label[[classVar]] <- as.factor(Train_label[[classVar]])
  fit <- rfsrc(formula = formula(paste0(classVar, "~.")),
               data = cbind(Train_set, Train_label[classVar]),
               ntree = 1000, nodesize = rf_nodesize,
               importance = T,
               proximity = T,
               forest = T)
  fit$subFeature = colnames(Train_set)
  return(fit)
}

RunGBM <- function(Train_set, Train_label, classVar){
  fit <- gbm(formula = Train_label[[classVar]] ~ .,
             data = as.data.frame(Train_set),
             distribution = 'bernoulli',
             n.trees = 10000,
             interaction.depth = 3,
             n.minobsinnode = 10,
             shrinkage = 0.001,
             cv.folds = 10,n.cores = 6)
  best <- which.min(fit$cv.error)
  fit <- gbm(formula = Train_label[[classVar]] ~ .,
             data = as.data.frame(Train_set),
             distribution = 'bernoulli',
             n.trees = best,
             interaction.depth = 3,
             n.minobsinnode = 10,
             shrinkage = 0.001, n.cores = 8)
  fit$subFeature = colnames(Train_set)
  return(fit)
}

RunXGBoost <- function(Train_set, Train_label, classVar){

  y <- Train_label[[classVar]]
  indexes = createFolds(y, k = 5, list = TRUE)
  CV <- unlist(lapply(indexes, function(pt){
    dtrain = xgb.DMatrix(
      data  = as.matrix(Train_set[-pt, , drop = FALSE]),
      label = y[-pt]
    )
    dtest = xgb.DMatrix(
      data  = as.matrix(Train_set[pt, , drop = FALSE]),
      label = y[pt]
    )
    watchlist <- list(train = dtrain, test = dtest)
    bst <- xgb.train(
      data     = dtrain,
      max.depth= 2,
      eta      = 1,
      nthread  = 2,
      nrounds  = 10,
      watchlist= watchlist,
      objective= "binary:logistic",
      verbose  = FALSE
    )
    which.min(bst$evaluation_log$test_logloss)
  }))
  nround <- as.numeric(names(which.max(table(CV))))
  fit <- xgboost(
    data     = as.matrix(Train_set),
    label    = y,
    max.depth= 2,
    eta      = 1,
    nthread  = 2,
    nrounds  = nround,
    objective= "binary:logistic",
    verbose  = FALSE
  )
  fit$subFeature = colnames(Train_set)
  return(fit)
}

RunNaiveBayes <- function(Train_set, Train_label, classVar){
  data <- cbind(Train_set, Train_label[classVar])
  data[[classVar]] <- as.factor(data[[classVar]])
  fit <- naiveBayes(eval(parse(text = paste(classVar, "~."))),
                    data = data)
  fit$subFeature = colnames(Train_set)
  return(fit)
}

quiet <- function(...) suppressMessages(eval(...))

ExtractVar <- function(fit){
  Feature <- quiet(switch(
    EXPR = class(fit)[1],
    "lognet" = rownames(coef(fit))[which(coef(fit)[, 1]!=0)],
    "glm" = names(coef(fit)),
    "svm.formula" = fit$subFeature,
    "train" = fit$coefnames,
    "glmboost" = names(coef(fit)[abs(coef(fit))>0]),
    "plsRglmmodel" = rownames(fit$Coeffs)[fit$Coeffs!=0],
    "rfsrc" = var.select(fit, verbose = F)$topvars,
    "gbm" = rownames(summary.gbm(fit, plotit = F))[summary.gbm(fit, plotit = F)$rel.inf>0],
    "xgb.Booster" = fit$subFeature,
    "naiveBayes" = fit$subFeature
  ))

  Feature <- setdiff(Feature, c("(Intercept)", "Intercept"))
  return(Feature)
}

CalPredictScore <- function(fit, new_data, type = "lp"){
  new_data <- new_data[, fit$subFeature]
  RS <- quiet(switch(
    EXPR = class(fit)[1],
    "lognet"      = predict(fit, type = 'response', as.matrix(new_data)),
    "glm"         = predict(fit, type = 'response', as.data.frame(new_data)),
    "svm.formula" = predict(fit, as.data.frame(new_data), probability = T),
    "train"       = predict(fit, new_data, type = "prob")[[2]],
    "glmboost"    = predict(fit, type = "response", as.data.frame(new_data)),
    "plsRglmmodel" = predict(fit, type = "response", as.data.frame(new_data)),
    "rfsrc"        = predict(fit, as.data.frame(new_data))$predicted[, "1"],
    "gbm"          = predict(fit, type = 'response', as.data.frame(new_data)),
    "xgb.Booster" = predict(fit, as.matrix(new_data)),
    "naiveBayes" = predict(object = fit, type = "raw", newdata = new_data)[, "1"]
  ))
  RS = as.numeric(as.vector(RS))
  names(RS) = rownames(new_data)
  return(RS)
}

PredictClass <- function(fit, new_data){
  new_data <- new_data[, fit$subFeature]
  label <- quiet(switch(
    EXPR = class(fit)[1],
    "lognet"      = predict(fit, type = 'class', as.matrix(new_data)),
    "glm"         = ifelse(test = predict(fit, type = 'response', as.data.frame(new_data))>0.5,
                           yes = "1", no = "0"),
    "svm.formula" = predict(fit, as.data.frame(new_data), decision.values = T),
    "train"       = predict(fit, new_data, type = "raw"),
    "glmboost"    = predict(fit, type = "class", as.data.frame(new_data)),
    "plsRglmmodel" = ifelse(test = predict(fit, type = 'response', as.data.frame(new_data))>0.5,
                            yes = "1", no = "0"),
    "rfsrc"        = predict(fit, as.data.frame(new_data))$class,
    "gbm"          = ifelse(test = predict(fit, type = 'response', as.data.frame(new_data))>0.5,
                            yes = "1", no = "0"),
    "xgb.Booster" = ifelse(test = predict(fit, as.matrix(new_data))>0.5,
                           yes = "1", no = "0"),
    "naiveBayes" = predict(object = fit, type = "class", newdata = new_data)
  ))
  label = as.character(as.vector(label))
  names(label) = rownames(new_data)
  return(label)
}

classVar = "outcome"
min.selected.var = 5

Variable = colnames(Train_set)
preTrain.method =  strsplit(methods, "\\+")
preTrain.method = lapply(preTrain.method, function(x) rev(x)[-1])
preTrain.method = unique(unlist(preTrain.method))

preTrain.var <- list()
set.seed(seed = 123)
for (method in preTrain.method){
  preTrain.var[[method]] = RunML(method = method,
                                 Train_set = Train_set,
                                 Train_label = Train_class,
                                 mode = "Variable",
                                 classVar = classVar)
}
preTrain.var[["simple"]] <- colnames(Train_set)

model <- list()
set.seed(seed = 123)
Train_set_bk = Train_set

for (method in methods){

  method_name = method
  method <- strsplit(method, "\\+")[[1]]

  if (length(method) == 1) method <- c("simple", method)

  Variable = preTrain.var[[method[1]]]
  Train_set = Train_set_bk[, Variable]
  Train_label = Train_class

  model[[method_name]] <- RunML(method = method[2],
                                Train_set = Train_set,
                                Train_label = Train_label,
                                mode = "Model",
                                classVar = classVar)

  if(length(ExtractVar(model[[method_name]])) <= min.selected.var) {
    model[[method_name]] <- NULL
  }
}

Train_set = Train_set_bk; rm(Train_set_bk)

  logisticmodel <- lapply(model, function(fit){

    tmp <- glm(formula = Train_class[[classVar]] ~ .,
               family = "binomial",
               data = as.data.frame(Train_set[, ExtractVar(fit)]))
    tmp$subFeature <- ExtractVar(fit)
    return(tmp)
  })


# Preserve the source evaluation order, using discovery data only.
methodsValid <- names(model)
RS_list <- list()
for (method in methodsValid) {
  RS_list[[method]] <- CalPredictScore(model[[method]], Train_set)
}
RS_mat <- as.data.frame(t(do.call(rbind, RS_list)))
write.table(RS_mat, file.path(OUT_DIR, "discovery_scores.tsv"),
            sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
Class_list <- list()
for (method in methodsValid) {
  Class_list[[method]] <- PredictClass(model[[method]], Train_set)
}
Class_mat <- as.data.frame(t(do.call(rbind, Class_list)))
write.table(Class_mat, file.path(OUT_DIR, "discovery_predictions.tsv"),
            sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
fea_list <- list()
for (method in methodsValid) {
  fea_list[[method]] <- ExtractVar(model[[method]])
}
fea_df <- lapply(model, function(fit) data.frame(ExtractVar(fit)))
fea_df <- do.call(rbind, fea_df)
fea_df$algorithm <- gsub("(.+)\\.(.+$)", "\\1", rownames(fea_df))
colnames(fea_df)[1] <- "features"
write.csv(fea_df, file.path(OUT_DIR, "algorithm_features.csv"), row.names = FALSE)
AUC_list <- list()
for (method in methodsValid) {
  score <- CalPredictScore(model[[method]], Train_set)
  AUC_list[[method]] <- if (length(unique(Train_class$outcome)) < 2) NA_real_ else
    as.numeric(auc(roc(Train_class$outcome, score)))
}
# discovery_all is the original pooled training AUC, not an out-of-fold estimate.
AUC_mat <- data.frame(Algorithm = names(AUC_list),
                     discovery_all = unlist(AUC_list), row.names = NULL)
AUC_mat <- AUC_mat[order(AUC_mat$discovery_all, decreasing = TRUE), ]
write.csv(AUC_mat, file.path(OUT_DIR, "algorithm_ranking.csv"), row.names = FALSE)
if (!nrow(AUC_mat) || all(is.na(AUC_mat$discovery_all)))
  stop("No evaluable discovery models remain.")
fea_sel <- fea_list[[AUC_mat$Algorithm[[1]]]]
write.table(fea_sel, file.path(OUT_DIR, "selected_features.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
saveRDS(model, file.path(OUT_DIR, "models.rds"))
saveRDS(logisticmodel, file.path(OUT_DIR, "logistic_models.rds"))
