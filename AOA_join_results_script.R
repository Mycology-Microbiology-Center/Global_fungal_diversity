#!/usr/bin/env Rscript

## Script to merge AOA chunks into a single raster

## Inputs:
# args[1] = Directrory with prediction chunks (with NAs, cellID and coordinates)
# args[2] = Input directory with *.RData chunks (results of `AOA_script.R`)
# args[3] = Name of the output file (without extension)

## NB. the file with prediction mask (`PredictionMask.tif`) should be in the current directory


args <- commandArgs(trailingOnly = TRUE)


library(terra)
library(data.table)
library(plyr)

library(doFuture)
registerDoFuture()
plan(multicore, workers = 10)          # will crash RStudio
options(future.globals.maxSize = 5e9)  # 5GB; default = 500 * 1024 ^ 2 = 500 MiB

OUTPUT <- args[3]

## List predicton chunks
PREDICTDIR <- args[1]
fls <- list.files(
  path = PREDICTDIR, pattern = ".RData", full.names = T, recursive = F)

cat("..Extracting AOAs from each chunk\n")

## Extract AOAs from each chunk
AOA <- alply(.data = fls, .margins = 1, .fun = function(x){
  # x <- fls[1]

  ## Load cell IDs from the predictors data
  xx <- readRDS(x)
  clz <- c("cellid", "LON", "LAT")
  xx <- xx[, ..clz]

  ## Load AOA
  aa_name <- file.path(
    args[2],
    gsub(pattern = ".RData", replacement = "_AOA.RData", x = basename(x))
    )
  aa <- readRDS(aa_name)

  res <- cbind(xx, data.table(DI = aa$DI, AOA = aa$AOA))

  attr(res, which = "threshold") <- attr(aa$AOA, "aoa_stats")$threshold
  return(res)

}, .progress = "text", .parallel = TRUE)

## Check thresholds
# ldply(.data = AOA, .fun = function(x){ data.frame(Tr = attr(x, which = "threshold")) })

threshold <- attr(AOA[[1]], which = "threshold")

cat("..Merging results into a single table\n")
AOA <- rbindlist(AOA)

cat("..Exporting results\n")
saveRDS(object = AOA,
  file = paste0(OUTPUT, ".RData"),
  compress = "xz")


### Re-create raster
cat("..Re-creating raster\n")

## Load prediction mask
cat("...Loading prediction mask\n")
RST <- rast("PredictionMask.tif")

## Drop raster vaues
cat("...Blanking the raster\n")
RST[ !is.na(RST) ] <- NA

## Insert values on a raster
cat("...Adding AOA values into the raster\n")
RST_AOA <- RST
RST_AOA[ AOA$cellid ] <- AOA$AOA
names(RST_AOA) <- "AOA"

cat("...Adding DI values into the raster\n")
RST_DI <- RST
RST_DI[ AOA$cellid ] <- AOA$DI
names(RST_DI) <- "DI"

cat("...Combining rasters\n")
RST_ALL <- c(RST_AOA, RST_DI)

## Export raster
cat("...Exporting raster\n")
writeRaster(
  x = RST_ALL,
  filename = paste0(OUTPUT, ".tif"),
  overwrite = TRUE,
  gdal=c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=8"))

cat("All done\n")
