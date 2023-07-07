#!/usr/bin/env Rscript

## Prepare data for dissimilarity index estimation
## Result = DI and var importances for AOA estimation

args <- commandArgs(trailingOnly = TRUE)

# args[1] = Input file with fitted model
# args[2] = Name of the output file


library(data.table)
library(CAST)
library(mlr3verse)
library(plyr)

## Load the trained mlr3 model
cat("..Loading a model\n")
MODELFILE <- args[1]         # e.g., "XGBoost_Models/GSMc_AllFungi_Model.RData"
MOD <- readRDS(MODELFILE)

## Extract feature importance
cat("..Extracting feature importance\n")
IMP <- data.table(
  Feature = names(MOD$learner$importance()),
  Importance = MOD$learner$importance())

## Match training and prediction data (+ remove Copernicus columns and EcM presence column)
cat("..Matching the data\n")
train_data <- attr(MOD, "tsk")$data()
train_data[, DependentVariable := NULL ]

cat("..Removing categorical variables\n")
if(any(grepl(pattern = "Copern_", x = colnames(train_data)))){
  clz <- grep(pattern = "Copern_", x = colnames(train_data), invert = T, value = T)
  clz <- grep(pattern = "EcM_Plant_Presence", x = clz, invert = T, value = T)
  train_data <- train_data[, ..clz]
  rm(clz)
}
clz <- colnames(train_data)


## Check if variables are present
if(any(! colnames(train_data) %in% IMP$Feature)){
  vrz <- paste(colnames(train_data)[ ! colnames(train_data) %in% IMP$Feature ], collapse = ", ")
  cat("...WARNING: some variables in train_data does not have importance values!\n")
  cat("... these variables will be excluded:\n")
  cat("   ", vrz, "\n")

  clz <- colnames(train_data)[ colnames(train_data) %in% IMP$Feature ]
  train_data <- train_data[, ..clz ]
}

## Subset feature importance to the predictors
IMP <- IMP[ Feature %in% clz ]

## Feature importances will be used as weights
imp_vars <- IMP$Feature
imp_wies <- t(IMP[, .(Importance)])
colnames(imp_wies) <- imp_vars

## Reorder columns (should be identical in all datasets and var weigths)
# setcolorder(preds_data, neworder = colnames(imp_wies))
setcolorder(train_data, neworder = colnames(imp_wies))

## Modified function for Dissimilarity Index of training data calculation
trainDI <- function(
  train = NULL,
  weight = NA,
  folds = NULL){

  train_backup <- train

  ## Scale training data
  train <- scale(train)

  ## Save scale param for later
  scaleparam <- attributes(train)

  # multiply train data with variable weights (from variable importance)
  train <- sapply(1:ncol(train),function(x){train[,x]*unlist(weight[x])})

  # make sure all variables have variance
  FAILS <- apply(train, 2, FUN=function(x){all(is.na(x))})
  if(any(FAILS)){
    cat("WARNING: There are some variables with NAs after normalization\n")
    FAILS <- colnames(train_backup)[ FAILS ]
  } else {
    FAILS <- NULL
  }

  # calculate average mean distance between training data
  trainDist_avrg <- c()
  trainDist_min <- c()

  for(i in seq(nrow(train))){

    # distance to all other training data (for average)
    trainDistAll <- FNN::knnx.dist(train, t(train[i,]), k = nrow(train))[-1]
    trainDist_avrg <- append(trainDist_avrg, mean(trainDistAll, na.rm = TRUE))

    # calculate distance to other training data:
    trainDist <- FNN::knnx.dist(t(matrix(train[i,])),train,k=1)
    trainDist[i] <- NA


    # mask of other folds
    if (!is.null(folds)){
      trainDist[folds==folds[i]] <- NA
    }

    trainDist_min <- append(trainDist_min, min(trainDist, na.rm = TRUE))
  }
  trainDist_avrgmean <- mean(trainDist_avrg, na.rm=TRUE)

  # Dissimilarity Index of training data -----
  TrainDI <- trainDist_min/trainDist_avrgmean

  # AOA Threshold ----
  thres <- grDevices::boxplot.stats(TrainDI)$stats[5]
  lower_thres <- grDevices::boxplot.stats(TrainDI)$stats[1]

  # Return: trainDI Object -------
  aoa_results = list(
    train = train_backup,
    weight = weight,
    variables = colnames(weight),
    catvars = NULL,
    scaleparam = scaleparam,
    trainDist_avrg = trainDist_avrg,
    trainDist_avrgmean = trainDist_avrgmean,
    trainDI = TrainDI,
    threshold = thres,
    lower_threshold = lower_thres
  )

  class(aoa_results) = "trainDI"
  attr(aoa_results, which = "FailedVariables") <- FAILS
  return(aoa_results)
}

setDF(train_data)

cat("..Preparing DI object\n")
trainDI_res <- trainDI(
  train = train_data,
  weight = imp_wies,
  folds = NULL)

## Export DI and var importances
cat("..Exporting results\n")
save(
  train_data, trainDI_res, imp_wies, imp_vars,
  file = args[2],
  compress = "xz")

