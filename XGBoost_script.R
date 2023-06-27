#!/usr/bin/Rscript

## Script to train XGBoost regression model (with autotuner)

## Usage:
# ./XGBoost_script.R \
#    --input "GSM_EcM.RData" \
#    --threads 32 \
#    --boruta 500 \
#    --autotuner 5000 \
#    --output "GSM_EcM_Model.RData"


## Input data format should be:
# column 1 = sample ID
# columns 2 and 3 = sample coordinates ("LAT", "LON")
# column 4 = dependent variable
# the other columns = predictors (all numeric)


############################################## Parse input parameters

## Check time
start_time <- Sys.time()

cat("Parsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse)) # don't say "Loading required package: optparse"

## Parse arguments
option_list <- list(
  make_option(c("-i", "--input"), action="store", default=NA, type='character', help="Input data"),
  make_option(c("-t", "--threads"), action="store", default=8L, type='integer', help="Number of CPUs to use"),
  make_option(c("-b", "--boruta"), action="store", default=500L, type='integer', help="Number of Boruta runs for feature pre-selection"),
  make_option(c("-a", "--autotuner"), action="store", default=10000L, type='integer', help="Number of autotuner evaluations"),
  make_option(c("-m", "--method"), action="store", default="random", type='character', help="Autotuner method (grid or random)"),
  make_option(c("-o", "--output"), action="store", default=NA, type='character', help="Output file name")
)
opt <- parse_args(OptionParser(option_list=option_list))


## Assign variables
INPUT <- opt$input
NCORES <- as.integer( opt$threads )
BORUTARUNS <- as.integer( opt$boruta )
METHOD <- opt$method
TUNEREVALUATIONS <- as.integer( opt$autotuner )
OUTPUT <- opt$output

## Log assigned variables
cat(paste("Input file: ", INPUT, "\n", sep=""))
cat(paste("Number of threads: ", NCORES, "\n", sep=""))
cat(paste("Number of Boruta runs: ", BORUTARUNS, "\n", sep=""))
cat(paste("Tuning method: ", METHOD, "\n", sep=""))
if(METHOD %in% "random"){
  cat(paste("Number of autotuner evaluations: ", TUNEREVALUATIONS, "\n", sep=""))
}
cat(paste("Output file: ", OUTPUT, "\n", sep=""))
cat("\n")



############################################## Load packages and data

cat("Loading packages...\n")

load_pckg <- function(pkg = "raster"){
  suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
  cat(paste(pkg, packageVersion(pkg), "\n"))
}

load_pckg("data.table")
load_pckg("mlr3verse")   # must be 0.13.3  remotes::install_github("mlr-org/mlr3@v0.13.3")
load_pckg("mlr3tuning")
load_pckg("mlr3learners")
load_pckg("mlr3misc")
load_pckg("paradox")     # to create parameter sets for hyperparameter tuning
load_pckg("Boruta")

cat("\n")

set.seed(14789)

lgr::get_logger("mlr3")$set_threshold("warn")

## Start future backend
# future::plan("multisession", workers = NCORES)

## Run both nested resampling loops (inner and outer) in parallel
future::plan(list(
  future::tweak("multisession", workers = 2),                  # outer loop (resampling)
  future::tweak("multisession", workers = round(NCORES/2))))   # inner loop (autotune)


## Load the data
cat("Loading input data...\n")
datt <- readRDS(INPUT)

## Remove rows with missing values
datt <- na.omit(datt)

############################################## Main workflow

## Rename dependent variable
DEPENDVAR <- colnames(datt)[4]
colnames(datt)[4] <- "DependentVariable"

idz <- "DependentVariable"
copern <- grep(pattern = "Copern_", x = colnames(datt), value = T)


## Run variable pre-selection using Boruta
## Do not take Copernicus variables into account
cat("Feature pre-selection with Boruta...\n")
feats <- Boruta(formula = DependentVariable ~ .,
  data = datt[, -which(colnames(datt) %in% c("Sample", "LAT", "LON", copern))],
  maxRuns = BORUTARUNS,
  num.threads = NCORES)

BORUTAFEATURES <- getSelectedAttributes(feats, withTentative = FALSE)
BORUTAIMPORTANCE <- attStats(feats)


## Subset the data
datt_subs <- datt[, unique(c(idz, copern, BORUTAFEATURES)) ]


cat("Setting up the model...\n")

## Create a regression task
tsk <- as_task_regr(datt_subs,
  target = "DependentVariable",
  id = "diversity")

## Specify the learner = XGboost
learner <- lrn("regr.xgboost")
set_threads(learner, n = 1)
#   learner$param_set
#   learner$param_set$values

## Define the inner CV (estimate model performance using 5-fold cross-validation)
resampling <- rsmp("cv", folds = 5)
measures <- msr("regr.rmse")


## Tune the model (AutoTuner automatically performs nested cross-validation)
# Termination Criterion = how long the tuning should run
# Tuner = which tuning method to use
# ParamSet = which space we might want to search through

# Define the search space
search_space = ps(
  eta = p_dbl(lower = 0.1, upper = 0.5),                       # learning rate (default, 0.3)
  max_depth = p_int(lower = 3L, upper = 10L),                  # depth of the tree (default, 6)
  min_child_weight = p_dbl(lower = 1, upper = 15),             # (default, 1)
  gamma = p_dbl(lower = 0, upper = 1),                         # regularization (default, 0 = no regul)
  subsample = p_dbl(lower = 0.7, upper = 1.0),                 # (default, 1)
  colsample_bytree = p_dbl( lower = 0.9, upper = 1),           # (default, 1)
  colsample_bylevel = p_dbl(lower = 0.5, upper = 1),           # (default, 1)
  nrounds = p_int(lower = 50L, upper = 400)                    # Number of trees in the ensemble  <NoDefault[3]> ?
  )

## Random search
if(METHOD %in% "random"){
  at <- auto_tuner(
    method = "random_search",
    learner = learner,
    resampling = resampling,          # resampling strategy = 5-fold cross-validation (inner loop)
    measure = measures,               # performance measure = RMSE
    search_space = search_space,
    term_evals = TUNEREVALUATIONS,    # number of random search iterations
    batch_size = NCORES)

  ## Train with AutoTuner
  cat("Training the model (with AutoTuner)...\n")
  at$train(tsk)
}

# ## Grid search
# if(METHOD %in% "grid"){
# 
#   grid <- generate_design_grid(
#     param_set = search_space,
#     resolution = 3)
# 
#   at <- auto_tuner(
#     method = "grid_search",
#     learner = learner,
#     resampling = resampling,          # resampling strategy = 5-fold cross-validation (inner loop)
#     measure = measures,               # performance measure = RMSE
#     param_set = grid,
#     batch_size = NCORES,
#     terminator = term("none")         # no terminator since all design points should be evaluated
#     )
# }


## Add metadata
attr(at, which = "DependentVariable") <- DEPENDVAR
attr(at, which = "BorutaSelectedFeatures") <- BORUTAFEATURES
attr(at, which = "BorutaFeatureImportance") <- BORUTAIMPORTANCE
attr(at, which = "tsk") <- tsk

## Export trained model
cat("Exporting results...\n")
saveRDS(object = at,
  file = OUTPUT,
  compress = "xz")


##############################################

## Check time
end_time <- Sys.time()

tmm <- as.numeric(difftime(end_time, start_time, units = "min"))
cat("\nElapsed time: ", tmm, " minutes\n")


cat("\n")
cat("Session info:\n")
sessionInfo()
cat("\n")
