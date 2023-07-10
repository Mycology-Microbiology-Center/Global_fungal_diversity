#!/usr/bin/env Rscript

## Script to estimate uncertainty for XGBoost predictions using cross-validation or bootstrap
## Step-1 - Predicting diversity (for individual chunks of data) using multiple models

## Input:
# `model` is a list with
# $task = TaskRegr
# $mod  = ResampleResult (with multiple models - e.g., MOD$mod$learners)

## Usage:
# ./XGBoost_prediction_bootstrap_script.R \
#    --model      "GSMc_EcM_Model_ForUncertainty.RData" \
#    --predchunk  "PredicitonChunks/Chunk_001.RData" \


############################################## Parse input parameters

## Check time
start_time <- Sys.time()

cat("Parsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse)) # don't say "Loading required package: optparse"

## Parse arguments
option_list <- list(
  make_option(c("-m", "--model"), action="store", default=NA, type='character', help="Trained model"),
  make_option(c("-p", "--predchunk"), action="store", default=NA, type='character', help="Data for prediction (single chunk)")
)
opt <- parse_args(OptionParser(option_list=option_list))


## Assign variables
MODEL      <- opt$model
PREDCHUNK  <- opt$predchunk

## Prepare output name
OUTPUT <- paste0(
  sub(pattern = "_Model_ForUncertainty.RData", replacement = "", x = basename(MODEL)),
  "__",
  basename(PREDCHUNK)
  )

## Log assigned variables
cat(paste("Model: ",               MODEL,     "\n", sep=""))
cat(paste("Data for prediction: ", PREDCHUNK, "\n", sep=""))
cat(paste("Output file name: ",    OUTPUT,    "\n", sep=""))
cat("\n")


############################################## Load packages and data

cat("Loading packages...\n")

load_pckg <- function(pkg = "raster"){
  suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
  cat(paste(pkg, packageVersion(pkg), "\n"))
}

load_pckg("data.table")
load_pckg("mlr3")
load_pckg("plyr")
# load_pckg("dplyr")

cat("\n")

set.seed(14789)

lgr::get_logger("mlr3")$set_threshold("warn")

## Load the model
cat("Loading the model...\n")
MOD <- readRDS(MODEL)

## Prediction function
my_predict <- function(chunkid = "PredicitonChunks/Chunk_001.RData"){

  cat("\nChunk ID: ", chunkid, "\n")

  ## Load the prediction chunk
  cat("..Loading data for prediction (single chunk)\n")
  PRD <- readRDS(chunkid)

  ## Fill missing Copenicus columns
  cat("..Fixing Copenicus codes\n")
  if(! "Copern_20"  %in% colnames(PRD)){ PRD[ , Copern_20  := 0 ] }
  if(! "Copern_30"  %in% colnames(PRD)){ PRD[ , Copern_30  := 0 ] }
  if(! "Copern_40"  %in% colnames(PRD)){ PRD[ , Copern_40  := 0 ] }
  if(! "Copern_50"  %in% colnames(PRD)){ PRD[ , Copern_50  := 0 ] }
  if(! "Copern_61"  %in% colnames(PRD)){ PRD[ , Copern_61  := 0 ] }
  if(! "Copern_62"  %in% colnames(PRD)){ PRD[ , Copern_62  := 0 ] }
  if(! "Copern_70"  %in% colnames(PRD)){ PRD[ , Copern_70  := 0 ] }
  if(! "Copern_90"  %in% colnames(PRD)){ PRD[ , Copern_90  := 0 ] }
  if(! "Copern_111" %in% colnames(PRD)){ PRD[ , Copern_111 := 0 ] }
  if(! "Copern_112" %in% colnames(PRD)){ PRD[ , Copern_112 := 0 ] }
  if(! "Copern_115" %in% colnames(PRD)){ PRD[ , Copern_115 := 0 ] }
  if(! "Copern_118" %in% colnames(PRD)){ PRD[ , Copern_118 := 0 ] }

  ## Prediction
  cat("..Prediction\n")

  ## Single model
  ## res <- MOD$mod$learners[[1]]$predict_newdata(newdata = PRD, task = MOD$task)

  ## Models from all bootrstrap iterations
  RES <- llply(
    .data = MOD$mod$learners,
    .fun = function(x){ 
      # x <- MOD$mod$learners[[1]]
      cat(".... One of the learners\n")
      pr <- x$predict_newdata(newdata = PRD, task = MOD$task)
      pr <- as.data.table(pr)
      pr[, c("row_ids", "truth") := NULL ]
      return(pr)
    })

  cat("..Processing results\n")

  ## Merge predictions into a single table
  RES <- dplyr::bind_cols(RES, .name_repair = "unique")
  setDT(RES)
  setnames(x = RES, old = colnames(RES), new = paste0("Pred_", 1:ncol(RES)))

  ## Estimate uncertainty
  cat("..Estimating uncertainty\n")
  RES[,  `:=`(
    Var = apply(.SD, 1, var, na.rm = T),
    IQR = apply(.SD, 1, IQR, na.rm = T)),
    .SDcols = patterns('Pred_') ]

  ## Prepare final result
  clz <- c("cellid", "LON", "LAT")
  RES <- cbind(PRD[, ..clz], RES)

  ## Add metadata
  setattr(RES, name = "PredChunk", value = chunkid)
  setattr(RES, name = "Model", value = MODEL)

  return(RES)
}


cat("Running prediction tasks...\n")
PREDS <- my_predict(PREDCHUNK)

cat("..Exporting predicions\n")
saveRDS(object = PREDS, file = OUTPUT, compress = "xz")

cat("All done!\n")


##############################################

## Check time
end_time <- Sys.time()

tmm <- as.numeric(difftime(end_time, start_time, units = "min"))
cat("\nElapsed time: ", tmm, " minutes\n")


cat("\n")
cat("Session info:\n")
sessionInfo()
cat("\n")
