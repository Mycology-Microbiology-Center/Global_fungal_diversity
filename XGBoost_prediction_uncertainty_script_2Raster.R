#!/usr/bin/Rscript

## Script to estimate uncertainty for XGBoost predictions using cross-validation or bootstrap
## Step-2 - aggreation of prediction prediction chunks into a single raster

## Usage:
# ./XGBoost_prediction_bootstrap_script.R \
#    --predictdir "Predicitons" \
#    --rastermask "PredictionMask.tif" \
#    --output     "GSMc_EcM_Model_Bootstrap"


############################################## Parse input parameters

## Check time
start_time <- Sys.time()

cat("Parsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse)) # don't say "Loading required package: optparse"

## Parse arguments
option_list <- list(
  make_option(c("-p", "--predictdir"), action="store", default=NA, type='character', help="Directory predicted data"),
  make_option(c("-r", "--rastermask"), action="store", default=NA, type='character', help="Raster template"),
  make_option(c("-o", "--output"),     action="store", default=NA, type='character', help="Output file preifx")
)
opt <- parse_args(OptionParser(option_list=option_list))


## Assign variables
PREDICTDIR <- opt$predictdir
RASTER     <- opt$rastermask
OUTPUT     <- opt$output

## Log assigned variables
cat(paste("Directory with prediction chunks: ",  PREDICTDIR, "\n", sep=""))
cat(paste("Raster template (prediction mask): ", RASTER,     "\n", sep=""))
cat(paste("Output file prefix: ",                OUTPUT,     "\n", sep=""))
cat("\n")


############################################## Load packages and data

cat("Loading packages...\n")

load_pckg <- function(pkg = "raster"){
  suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
  cat(paste(pkg, packageVersion(pkg), "\n"))
}

load_pckg("data.table")
load_pckg("terra")
load_pckg("plyr")

cat("\n")

set.seed(14789)


## Load prediction mask
cat("Loading prediction mask...\n")
RST <- rast(RASTER)

cat("Loading predictions...\n")
## Files with predictions
fls <- list.files(path = PREDICTDIR, pattern = ".RData", full.names = T, recursive = F)
cat("..Files detected: ", length(fls), "\n")
PREDS <- alply(.data = fls,
  .margins = 1,
  .fun = readRDS,
  .progress = "text")

cat("Aggregating prediction chunks\n")
PREDS <- rbindlist(PREDS)

## Convert variance into SD
PREDS[ , SD := sqrt(Var) ]

cat("Exporting aggregated predicions\n")
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
