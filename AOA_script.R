#!/usr/bin/env Rscript

## Script to estimate area of applicability (AOA)

## Usage:
# ./AOA_script.R \
#    --grid "PredicitonChunks/Chunk_001.RData" \
#    --data "trainDI__GSMc_All.RData" \
#    --threads 10


############################################## Parse input parameters

## Check time
start_time <- Sys.time()

cat("Parsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse)) # don't say "Loading required package: optparse"

## Parse arguments
option_list <- list(
  make_option(c("-g", "--grid"), action="store", default=NA, type='character', help="Grid for prediction (data.table)"),
  make_option(c("-d", "--data"), action="store", default=NA, type='character', help="Training data, DI, and variable importances"),
  make_option(c("-t", "--threads"), action="store", default=10L, type='integer', help="Number of CPU threads")
  )
opt <- parse_args(OptionParser(option_list=option_list))

## Validation
# if(is.na(opt$input)){
#   cat("Input file is not specified.\n", file=stderr())
#   stop()
# }

## Assign variables
GRID <- opt$grid
DATA <- opt$data
THREADS <- opt$threads

## Prepare output file name
OUTPUT <- gsub(pattern = ".RData", replacement = "_AOA.RData", x = basename(GRID))

## Log assigned variables
cat(paste("Grid data: ", GRID, "\n", sep=""))
cat(paste("Training data: ", DATA, "\n", sep=""))
cat(paste("CPU threads: ", THREADS, "\n", sep=""))
cat(paste("Output file: ", OUTPUT, "\n", sep=""))

cat("\n")

############################################## Load packages

cat("Loading packages...\n")

load_pckg <- function(pkg = "raster"){
  suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
  cat(paste(pkg, packageVersion(pkg), "\n"))  
}

load_pckg("data.table")
load_pckg("CAST")
load_pckg("doParallel")
load_pckg("parallel")

## Set up the cluster
cat("Setting up the cluster for multi-threading..\n")
CLUSTER <- makeCluster(THREADS, type="FORK")
registerDoParallel(CLUSTER)


cat("\n")



############################################## Main workflow

## Load prediction chunk without NA values
## Variables were subsetted
cat("Loading grid data..\n")
preds_data <- readRDS(GRID)

## Load training data, DI, and variables importances
cat("Loading training data..\n")
load(DATA)
#   "imp_vars", "imp_wies", "train_data", "trainDI_res"

## Data validation
if(any(!colnames(train_data) %in% colnames(preds_data))){
  cat("WARNING: some of training data variables are missing from the predictors!\n")
  miss <- paste(colnames(train_data)[ !colnames(train_data) %in% colnames(preds_data)], collapse = ", ")
  cat(miss, "\n")
}

if(any(!colnames(preds_data) %in% colnames(train_data))){
  cat("Some of the predictors are not in the training data (that's OK).\n")
  miss2 <- paste(colnames(preds_data)[ !colnames(preds_data) %in% colnames(train_data)], collapse = ", ")
  cat(miss2, "\n")
}

## Re-order variables
cat("Re-ordering and subsetting predictors..\n")
clz <- colnames(train_data)
preds_data <- preds_data[, ..clz]


## Estimate the Area of Applicability (AOA)
# https://github.com/HannaMeyer/CAST/blob/master/R/aoa.R
cat("Estimating AOA..\n")
AOA <- aoa(newdata = preds_data,
           model = NA,
           trainDI = trainDI_res,
           train = train_data,                # as.data.frame(task$data()),
           variables = imp_vars,              # task$feature_names,
           weight = imp_wies,                 # data.frame(t(lrn$importance())),
           folds = NULL,                      # folds = rsmp_cv$instance[order(row_id)]$fold
           cl = CLUSTER)

cat("Stopping the cluster..\n")
stopCluster(CLUSTER)

cat("Saving results..\n")
saveRDS(object = AOA, file = OUTPUT, compress = "xz")

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

