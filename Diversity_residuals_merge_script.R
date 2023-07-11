#!/usr/bin/env Rscript

## Script to merge chunks of diversity residuals
## and re-create a full raster

## Usage:
# ./Diversity_residuals_merge_script.R \
#    --resids "./GSMc_AgarNM/" \
#    --model "GSMc_AgarNM__Raster.tif" \
#    --predictionmask "PredictionMask.tif" \
#    --output "GSMc_AgarNM"


############################################## Parse input parameters

## Check time
start_time <- Sys.time()

cat("Parsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse))

## Parse arguments
option_list <- list(
  make_option(c("-r", "--resids"), action="store", default=NA, type='character', help="Path to the residuals chunks"),
  make_option(c("-m", "--model"), action="store", default=NA, type='character', help="Model-based predictions (tif)"),
  make_option(c("-p", "--predictionmask"), action="store", default=NA, type='character', help="Prediction mask (tif)"),
  make_option(c("-o", "--output"), action="store", default="Out", type='character', help="Output file prefix")
)
opt <- parse_args(OptionParser(option_list=option_list))


## Assign variables
RESIDSPATH <- opt$resids
MODEL <- opt$model
MASK <- opt$predictionmask
OUTPUT <- opt$output

## Log assigned variables
cat(paste("Path to the residuals chunks: ", RESIDSPATH, "\n", sep=""))
cat(paste("Model-based predictions: ", MODEL, "\n", sep=""))
cat(paste("Prediction mask: ", MASK, "\n", sep=""))
cat(paste("Output file prefix: ", OUTPUT, "\n", sep=""))

cat("\n")

############################################## Load packages

cat("Loading packages...\n")

load_pckg <- function(pkg = "raster"){
  suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
  cat(paste(pkg, packageVersion(pkg), "\n"))  
}

load_pckg("data.table")
load_pckg("plyr")
load_pckg("terra")

cat("\n")



############################################## Main workflow


## Load and merge interpolated residuals
cat("Loading residuals..\n")
fls <- list.files(path = RESIDSPATH, pattern = ".RData", full.names = TRUE, recursive = FALSE)
RESIDS <- alply(.data = fls, .margins = 1, .fun = readRDS, .progress = "text")
RESIDS <- rbindlist(RESIDS)

## Load prediction mask
cat("Loading prediction mask..\n")
RST <- rast(MASK)
RST[ !is.na(RST) ] <- NA                    # Drop raster vaues
RST[ RESIDS$cellid ] <- RESIDS$Residuals    # Insert predicted values on a raster

## Add variable name to the raster
names(RST) <- "Residuals"

## Export raster
cat("Exporting interpolated residuals raster..\n")
writeRaster(
  x = RST,
  filename = paste0(OUTPUT, "__2.Residuals.tif"),
  overwrite=TRUE,
  gdal=c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=8"))

## Load model-based predictions
cat("Loading model-based predictions..\n")
MRS <- rast(MODEL)

## Add residuals to the model (a.k.a., regression-krigging approach)
RKG <- MRS + RST

## Export raster
cat("Exporting final raster...\n")
writeRaster(
  x = RKG,
  filename = paste0(OUTPUT, "__3.RegrKrig.tif"),
  overwrite=TRUE,
  gdal=c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=8"))


###### Plots

colors <- rev( colorspace::divergingx_hcl(n = 40, "RdYlBu") )

## Plot residuals
cat("Plotting residuals...\n")

jpeg(filename = paste0(OUTPUT, "__2.Residuals.jpg"), width = 20000, height = 10000, units = "px")
  par(oma = c(2,1,0,0) + 0.1, mar = c(0,0,1,1) + 0.1)

  plot(RST, col = colors, maxcell = 10000000)

dev.off()


## Plot final raster
cat("Plotting final raster...\n")

jpeg(filename = paste0(OUTPUT, "__3.RegrKrig.jpg"), width = 20000, height = 10000, units = "px")
  par(oma = c(2,1,0,0) + 0.1, mar = c(0,0,1,1) + 0.1)

  plot(RKG, col = colors, maxcell = 10000000)

dev.off()


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
