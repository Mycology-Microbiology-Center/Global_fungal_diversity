#!/usr/bin/env Rscript

## Script to estimate Shapley values for XGBoost models
## "Explanation Level Uncertainty of Sequential Variable Attribution"

## Usage:
# ./XGBoost_Shapley_script.R \
#    --input "ShapInp_2.RData" \
#    --explainer "GSMc_Pathog.RData" \
#    --threads 8

## `--threads` arg is not currently used


############################################## Parse input parameters

## Check time
start_time <- Sys.time()

cat("Parsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse)) # don't say "Loading required package: optparse"

## Parse arguments
option_list <- list(
  make_option(c("-i", "--input"), action="store", default=NA, type='character', help="Input data (data.table with predictors)"),
  make_option(c("-e", "--explainer"), action="store", default=NA, type='character', help="DALEX explainer"),
  make_option(c("-t", "--threads"), action="store", default=8L, type='integer', help="Number of CPUs to use")

)
opt <- parse_args(OptionParser(option_list=option_list))


## Assign variables
INPUT <- opt$input
EXPLAINER <- opt$explainer
NCORES <- as.integer( opt$threads )

OUTPUT <- sub(pattern = "ShapInp_",
              replacement = "ShapOUT_",
              x = basename(INPUT))

## Log assigned variables
cat(paste("Input file: ",  INPUT,     "\n", sep=""))
cat(paste("Explainer: ",   EXPLAINER, "\n", sep=""))
cat(paste("Output file: ", OUTPUT,    "\n", sep=""))
cat(paste("Number of threads: ", NCORES, "\n", sep=""))
cat("\n")


############################################## Load packages and data

cat("Loading packages...\n")

load_pckg <- function(pkg = "raster"){
  suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
  cat(paste(pkg, packageVersion(pkg), "\n"))
}

load_pckg("data.table")
load_pckg("xgboost")      # 1.6.0.1
load_pckg("fastshap")     # 0.0.7           remotes::install_github("bgreenwell/fastshap")
load_pckg("plyr")
# load_pckg("mlr3verse")  # must be 0.13.3  remotes::install_github("mlr-org/mlr3@v0.13.3")
# load_pckg("DALEX")      # 2.4.2           install.packages("DALEX")
# load_pckg("DALEXtra")   # 2.2.1           install.packages("DALEXtra")
# load_pckg("iBreakDown") # 2.0.1

cat("\n")

set.seed(14789)

## Disable mlr3 parallelization
# options(
#   parallelMap.default.mode        = "local",
#   parallelMap.default.cpus        = 1,
#   parallelMap.default.show.info   = FALSE
# )

## Set up a local cluster
# if(NCORES > 1){
#   cat("..Setting up a local cluster with ", NCORES, "cores\n")
#   library(doFuture)
#   registerDoFuture()
#   plan(multicore, workers = NCORES)       # will crash RStudio
#   options(future.globals.maxSize = 30e9)  # 30GB
#   parall <- TRUE
# } else {
#   parall <- FALSE
# }


## Load explainer
cat("Loading DALEX explainer...\n")
expl <- readRDS(EXPLAINER)
expl$PROFILES <- NULL

## Load environmental data
cat("Loading environmental data...\n")
ENV <- readRDS(INPUT)


cat("Checking variables...\n")
predz <- colnames(expl$EXPLAINER$data)
missing_predz <- predz[ ! predz %in% colnames(ENV) ]
if(length(missing_predz) > 0){
  ENV[, c(missing_predz) := NA ]
}


############################################## Main workflow

cat("Running analysis...\n")

## Subset data
## Matrix must containing ONLY the feature columns from the training data!
cat("..Preparing data\n")
setDF(ENV)
xx <- ENV[ , expl$EXPLAINER$model$model$feature_names ]

## Estimate Shapley values with `fastshap`
## Very fast, but approximate (nevertheless, values are quite similar to `iBreakDown`)
cat("..Estimating Shapley values\n")
RES <- fastshap::explain(
  object = expl$EXPLAINER$model$model,
  X = as.matrix(xx),
  pred_wrapper = pfun, exact = TRUE)

cat("..Preparing results\n")

RES <- as.data.frame(RES)

## Merge predictors and shap values
setDT(RES)
setDT(ENV)

colnames(ENV)[ ! colnames(ENV) %in% "cellid" ] <- paste0(
  "Env__", colnames(ENV)[ ! colnames(ENV) %in% "cellid" ])

colnames(RES) <- paste0("Shap__", colnames(RES))

RES <- cbind(ENV, RES)


## Export results
cat("Exporting results...\n")
saveRDS(object = RES,
  file = OUTPUT,
  compress = "xz")

cat("\nAll done\n")

##############################################

## Check time
end_time <- Sys.time()

tmm <- as.numeric(difftime(end_time, start_time, units = "min"))
cat("\nElapsed time: ", tmm, " minutes\n")

cat("\n")
cat("Session info:\n")
sessionInfo()
cat("\n")

