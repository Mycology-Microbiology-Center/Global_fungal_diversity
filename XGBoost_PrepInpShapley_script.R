#!/usr/bin/Rscript

## Script to prepare input data for Shaply values estimation
## chunks (10 000 regularly-spaced points per ecoregion)

## Inputs:
# args[1] = Ecoregion IDs (range, e.g., "1:10")

## Inputs:
# - Ecoregion polygons in sf format (`Ecoregions2017.RData`)
# - PredictionMask.tif
# - Directory with prediction chunks in Parquet format (`PredicitonChunks_parquet2`)

## The number of CPU threads is hardcoded to 10.


args <- commandArgs(trailingOnly = TRUE)


load_pckg <- function(pkg = "raster"){
  suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
  cat(paste(pkg, packageVersion(pkg), "\n"))
}

load_pckg("terra")
load_pckg("sf")
load_pckg("data.table")
load_pckg("arrow")
load_pckg("dplyr")

cat("Chunk: ", args[1], "\n")

## Split the range
spl <- strsplit(x = args[1], split = ":")[[1]]
MINID <- as.integer(spl[1])
MAXID <- as.integer(spl[2])
cat("Min eco ID: ", MINID, "\n")
cat("Max eco ID: ", MAXID, "\n")

set.seed(14789)

set_cpu_count(10)           # for libarrow
setDTthreads(threads = 10)  # for data.table

## Load reference raster (for cellids)
RST <- rast("PredictionMask.tif")

## Load ecoregion shapefiles
ECO17 <- readRDS("Ecoregions2017.RData")

## Load environmental data (arrow, WITH `cellid` column)
ENV <- open_dataset("PredicitonChunks_parquet/")

## Main function
prep_shap_inp <- function(
  ecoID = 1,
  npoints = 10000){

  ## Select ecoregion
  eco <- ECO17[ ecoID, ]

  ## Sample N points
  ppt <- suppressMessages(
    st_sample(eco, size = npoints, type = "regular")
    )
  st_crs(ppt) <- st_crs(eco)

  ## Point coordinates
  ppt <- st_coordinates(ppt)

  ## Get cell IDs from the raster
  ppt <- data.table(
    LON = ppt[,1],
    LAT = ppt[,2],
    cellid = cellFromXY(object = RST, xy = ppt))

  ppt <- unique(ppt, by = "cellid")

  ## Extract predictors by cellid (arrow)
  res <- ENV %>% filter(cellid %in% ppt$cellid) %>%  collect()
  setDT(res)

  if(nrow(res) > 0){

    ## Replace missing values
    res[, as.data.table(apply(.SD, 2, function(x){ x[is.na(x)] <- mean(x, na.rm=T); x })) ]

    ## Add ecoregion attributes
    res[, Olson_ECO_NAME   := eco$Olson_ECO_NAME   ]
    res[, Olson_BIOME_NAME := eco$Olson_BIOME_NAME ]
    res[, Olson_REALM      := eco$Olson_REALM      ]

    ## Export results
    saveRDS(object = res,
      file = file.path("Input", paste0("ShapInp_", ecoID, ".RData")),
      compress = "xz")

    ## Clean up
    rm(ppt, res, eco)
  } else {
    cat("\nWARNING: no points extracted for ", ecoID, " - ", eco$Olson_ECO_NAME, "\n")
    rm(ppt, eco)
  }
  gc()
}

plyr::a_ply(
  .data = MINID:MAXID,
  .margins = 1, .fun = function(x){
  try( prep_shap_inp(ecoID = x) )
  }, .progress = "text", .parallel = FALSE)

