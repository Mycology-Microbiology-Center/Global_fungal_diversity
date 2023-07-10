#!/usr/bin/Rscript

## Script to estimate uncertainty for XGBoost predictions using cross-validation or bootstrap

## Input:
# `model` is a list with
# $task = TaskRegr
# $mod  = ResampleResult (with multiple models - e.g., MOD$mod$learners)

## Usage:
# ./XGBoost_prediction_bootstrap_script.R \
#    --model "GSM_EcM_Model_ForUncertainty.RData" \
#    --predictdir "PredicitonChunks" \
#    --rastermask "PredictionMask.tif" \
#    --threads 10 \
#    --output "GSM_EcM_Model_Bootstrap"


############################################## Parse input parameters

## Check time
start_time <- Sys.time()

cat("Parsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse)) # don't say "Loading required package: optparse"

## Parse arguments
option_list <- list(
  make_option(c("-m", "--model"), action="store", default=NA, type='character', help="Trained model"),
  make_option(c("-p", "--predictdir"), action="store", default=NA, type='character', help="Directory with data for predictions"),
  make_option(c("-r", "--rastermask"), action="store", default=NA, type='character', help="Raster template"),
  make_option(c("-t", "--threads"), action="store", default=10L, type='integer', help="Number of CPUs to use"),
  make_option(c("-o", "--output"), action="store", default=NA, type='character', help="Output file preifx")
)
opt <- parse_args(OptionParser(option_list=option_list))


## Assign variables
MODEL      <- opt$model
PREDICTDIR <- opt$predictdir
RASTER     <- opt$rastermask
NCORES     <- as.integer( opt$threads )
OUTPUT     <- opt$output

## Log assigned variables
cat(paste("Model: ", MODEL, "\n", sep=""))
cat(paste("Directory with prediction chunks: ", PREDICTDIR, "\n", sep=""))
cat(paste("Raster template (prediction mask): ", RASTER, "\n", sep=""))
cat(paste("Number of threads: ", NCORES, "\n", sep=""))
cat(paste("Output file prefix: ", OUTPUT, "\n", sep=""))
cat("\n")


############################################## Load packages and data

cat("Loading packages...\n")

load_pckg <- function(pkg = "raster"){
  suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
  cat(paste(pkg, packageVersion(pkg), "\n"))
}

load_pckg("data.table")
load_pckg("mlr3")
load_pckg("terra")
load_pckg("plyr")
# load_pckg("dplyr")

cat("\n")

set.seed(14789)

lgr::get_logger("mlr3")$set_threshold("warn")

## Start future backend
if(NCORES > 1){
  cat("Starting future backend...\n")

  load_pckg("doFuture")
  registerDoFuture()
  # plan(multisession, workers = NCORES)    # for RStudio
  plan(multicore, workers = NCORES)         # will crash RStudio
  options(future.globals.maxSize = 5e9)     # 5GB; default = 500 * 1024 ^ 2 = 500 MiB
  options(doFuture.rng.onMisuse = "ignore")

  parall <- TRUE
} else {
  parall <- FALSE
}


## Load the model
cat("Loading the model...\n")
MOD <- readRDS(MODEL)

## Files with data for prediction
fls <- list.files(path = PREDICTDIR, pattern = ".RData", full.names = T, recursive = F)

## Load prediction mask
cat("Loading prediction mask...\n")
RST <- rast(RASTER)


############################################## Main workflow

## Prediction function
my_predict <- function(chunkid = "PredicitonChunks/Chunk_001.RData"){

  cat("\n..Chunk ID: ", chunkid, "\n")

  ## Load the prediction chunk
  cat("...Loading prediction chunk\n")
  # PRD <- readRDS("PredicitonChunks/Chunk_001.RData")
  PRD <- readRDS(chunkid)

  ## Fill missing Copenicus columns
  cat("...Fixing Copenicus codes\n")
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
  cat("...Prediction\n")

  ## Models from all bootrstrap iterations
  RES <- llply(
    .data = MOD$mod$learners,
    .fun = function(x){ 
      # x <- MOD$mod$learners[[1]]
      pr <- x$predict_newdata(newdata = PRD, task = MOD$task)
      pr <- as.data.table(pr)
      pr[, c("row_ids", "truth") := NULL ]
      return(pr)
    })

  cat("...Processing results\n")

  ## Merge predictions into a single table
  # RES <- do.call(cbind, RES)
  RES <- dplyr::bind_cols(RES, .name_repair = "unique")
  setDT(RES)
  setnames(x = RES, old = colnames(RES), new = paste0("Pred_", 1:ncol(RES)))

  ## Estimate uncertainty
  RES[,  `:=`(
    Var = apply(.SD, 1, var, na.rm = T),
    IQR = apply(.SD, 1, IQR, na.rm = T)),
    .SDcols = patterns('Pred_') ]

  ## Prepare final result
  clz <- c("cellid", "LON", "LAT")
  RES <- cbind(PRD[, ..clz], RES)

  ## Add metadata
  setattr(RES, name = "DependentVariable", value = attr(MOD, which = "DependentVariable"))
  setattr(RES, name = "PredChunk", value = chunkid)
  setattr(RES, name = "Model", value = MODEL)

  return(RES)
}


cat("Running prediction tasks...\n")
PREDS <- alply(.data = fls,
  .margins = 1,
  .fun = my_predict,
  .progress = "text",
  .parallel = parall)


cat("..Aggregating prediction chunks\n")
PREDS <- rbindlist(PREDS)

## Convert variance into SD
PREDS[ , SD := sqrt(Var) ]

cat("..Exporting predicions\n")
saveRDS(object = PREDS,
  file = paste0(OUTPUT, "__UnsertaintyPredictions.RData"),
  compress = "xz")


## Re-create raster
cat("Re-creating the raster...\n")

## Drop raster vaues
cat("..Removing raster values\n")
RST[ !is.na(RST) ] <- NA

## Insert predicted values on a raster
cat("..Inserting predicted values on a raster\n")

RSS <- copy(RST)
RSS[ PREDS$cellid ] <- PREDS$SD
names(RSS) <- "SD"

RSI <- copy(RST)
RSI[ PREDS$cellid ] <- PREDS$IQR
names(RSI) <- "IQR"

## Merge rasters
RST <- c(RSS, RSI)

## Export raster
cat("..Exporting the raster\n")

writeRaster(x = RST,
  filename = paste0(OUTPUT, "__UncertaintyRaster.tif"),
  overwrite = TRUE,
  gdal=c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=8"))

cat("..Plotting the raster\n")

# colors <- rev( scico::scico(40, palette = 'nuuk') )
colors <- rev(
  c("#04598C", "#0E5B89", "#175D87", "#1D5F85", "#256083", "#2B6380", 
    "#336780", "#3A6A80", "#436E81", "#4B7281", "#557884", "#5E7C86", 
    "#67828A", "#71878D", "#7A8D91", "#839394", "#8B9795", "#949D98", 
    "#9AA199", "#A1A599", "#A6A998", "#ABAD97", "#AFAF96", "#B3B393", 
    "#B5B591", "#B8B88D", "#BABA8B", "#BEBD88", "#C0C085", "#C3C282", 
    "#C7C681", "#CAC980", "#CFCF7F", "#D5D481", "#DCDC86", "#E3E38C", 
    "#EBEB95", "#F2F29E", "#F8F8A8", "#FEFEB2"))

## JPG
cat("... SD\n")
jpeg(filename = paste0(OUTPUT, "__Uncertainty_SD.jpg"), width = 20000, height = 10000, units = "px")
  par(oma = c(2,1,0,0) + 0.1, mar = c(0,0,1,1) + 0.1)

  plot(RST[[1]], col = colors, maxcell = 10000000)

dev.off()

cat("... IQR\n")
jpeg(filename = paste0(OUTPUT, "__Uncertainty_IQR.jpg"), width = 20000, height = 10000, units = "px")
  par(oma = c(2,1,0,0) + 0.1, mar = c(0,0,1,1) + 0.1)

  plot(RST[[2]], col = colors, maxcell = 10000000)

dev.off()

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
