#!/usr/bin/Rscript

## Script to apply XGBoost regression model on the environmental data 
## to predict the expected diversity,
## estimat diversity residuals,
## and create and plot the raster with expected diversity

## Usage:
# ./XGBoost_prediction_script.R \
#    --model "GSM_EcM_Model.RData" \
#    --input "Pipe_inp/GSM_EcM.RData" \
#    --predictdir "PredicitonChunks" \
#    --rastermask "PredictionMask.tif" \
#    --threads 10 \
#    --output "GSM_EcM_Model"


############################################## Parse input parameters

## Check time
start_time <- Sys.time()

cat("Parsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse)) # don't say "Loading required package: optparse"

## Parse arguments
option_list <- list(
  make_option(c("-m", "--model"), action="store", default=NA, type='character', help="Trained model"),
  make_option(c("-i", "--input"), action="store", default=NA, type='character', help="Input data"),
  make_option(c("-p", "--predictdir"), action="store", default=NA, type='character', help="Directory with data for predictions"),
  make_option(c("-r", "--rastermask"), action="store", default=NA, type='character', help="Raster template"),
  make_option(c("-t", "--threads"), action="store", default=10L, type='integer', help="Number of CPUs to use"),
  make_option(c("-o", "--output"), action="store", default=NA, type='character', help="Output file preifx")
)
opt <- parse_args(OptionParser(option_list=option_list))


## Assign variables
MODEL <- opt$model
INPUT <- opt$input
PREDICTDIR <- opt$predictdir
RASTER <- opt$rastermask
NCORES <- as.integer( opt$threads )
OUTPUT <- opt$output

## Log assigned variables
cat(paste("Trained model: ", MODEL, "\n", sep=""))
cat(paste("Input data: ", INPUT, "\n", sep=""))
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
load_pckg("mlr3verse")
load_pckg("terra")
load_pckg("plyr")

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

  parall <- TRUE
} else {
  parall <- FALSE
}


## Load the trained model
cat("Loading the trained model...\n")
MOD <- readRDS(MODEL)

## File with source data (for residuals)
cat("Loading the input data...\n")
SRC <- readRDS(INPUT)
setDT(SRC)

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
  RES <- MOD$predict_newdata(newdata = PRD)
  
  cat("...Processing results\n")
  RES <- as.data.table(RES)

  ## Prepare final result
  clz <- c("cellid", "LON", "LAT")
  REZ <- PRD[, ..clz]
  REZ$Prediction <- RES$response

  ## Add metadata
  setattr(REZ, name = "DependentVariable", value = attr(MOD, which = "DependentVariable"))
  setattr(REZ, name = "PredChunk", value = chunkid)
  setattr(REZ, name = "Model", value = MODEL)

  return(REZ)
}

cat("Running prediction tasks...\n")
PREDZ <- alply(.data = fls, .margins = 1, .fun = my_predict, .progress = "text", .parallel = parall)

## Extract the name of dependent variable
DependentVariable <- attr(PREDZ[[1]], which = "DependentVariable")

cat("..Aggregating prediction chunks\n")
PREDS <- rbindlist(PREDZ)

cat("..Exporting predicions\n")
saveRDS(object = PREDS, file = paste0(OUTPUT, "__Predictions.RData"), compress = "xz")


## Get residuals
cat("Estimating residuals...\n")
tmp <- MOD$predict_newdata(newdata = SRC)
tmp <- as.data.table(tmp)
RESIDUALS <- SRC[, 1:4]
colnames(RESIDUALS)[4] <- "Observed"
RESIDUALS$Predicted <- tmp$response
setattr(x = RESIDUALS, name = "DependentVariable", value = colnames(SRC)[4])

RESIDUALS[, Residuals := Observed - Predicted ]

cat("..Exporting residuals\n")
saveRDS(object = RESIDUALS, file = paste0(OUTPUT, "__Residuals.RData"), compress = "xz")


## Re-create raster
cat("Re-creating the raster...\n")

## Drop raster vaues
cat("..Removing raster values\n")
RST[ !is.na(RST) ] <- NA

## Insert predicted values on a raster
cat("..Inserting predicted values on a raster\n")
RST[ PREDS$cellid ] <- PREDS$Prediction

## Add variable name to the raster
names(RST) <- DependentVariable

## Export raster
cat("..Exporting the raster\n")

writeRaster(x = RST,
  filename = paste0(OUTPUT, "__Raster.tif"),
  overwrite = TRUE,
  gdal=c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=8"))

cat("..Plotting the raster\n")

colors <- rev( colorspace::divergingx_hcl(n = 40, "RdYlBu") )

## JPG
jpeg(filename = paste0(OUTPUT, "__Raster.jpg"), width = 20000, height = 10000, units = "px")
  par(oma = c(2,1,0,0) + 0.1, mar = c(0,0,1,1) + 0.1)

  plot(RST, col = colors, maxcell = 10000000)

dev.off()

## PDF
# pdf(file = paste0("EcM_map", ".pdf"), width = 30, height = 15, useDingbats = FALSE)
#       par(oma = c(2,1,0,0) + 0.1, mar = c(0,0,1,1) + 0.1)
# 
#   plot(RST, col = colors, maxcell = 5000000)
# 
# dev.off()


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
