#!/usr/bin/env Rscript

## Estimate diversity residuals on a grid
## Input:
## - model-based residuals (data.table with columns  "Sample", "LAT", "LON", "Observed", "Predicted", "Residuals")
## - grid points for prediction (data.table with columns "cellid", "LAT", "LON")


## Usage:
# ./Diversity_residuals_script.R \
#    --data "Residuals_ECM.RData" \
#    --grid "PredicitonChunks/Chunk_001.RData" \
#    --idwpower 1.4 \
#    --output Results.RData


# IDW: Higher power value (`idwpower`) puts more emphasis on the nearest points.
#      Thus, nearby data will have the most influence, and the surface will have more detail (be less smooth).
#      As the power increases, the interpolated values begin to approach the value of the nearest sample point.
#      Specifying a lower value for power will give more influence to surrounding points that are farther away,
#      resulting in a smoother surface.

############################################## Parse input parameters

## Check time
start_time <- Sys.time()

cat("Parsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse)) # don't say "Loading required package: optparse"

## Parse arguments
option_list <- list(
  make_option(c("-d", "--data"), action="store", default=NA, type='character', help="File with observed residuals (data.table)"),
  make_option(c("-g", "--grid"), action="store", default=NA, type='character', help="Grid for prediction (data.table)"),
  make_option(c("-p", "--idwpower"), action="store", default=2, type='double', help="Inverse distance weighting power"),
  make_option(c("-o", "--output"), action="store", default="Out.RData", type='character', help="Output file name")
)
opt <- parse_args(OptionParser(option_list=option_list))

## Validation
# if(is.na(opt$input)){
#   cat("Input file is not specified.\n", file=stderr())
#   stop()
# }

## Assign variables
DATA <- opt$data
GRID <- opt$grid
IDWPOWER <- opt$idwpower
OUTPUT <- opt$output

## Log assigned variables
cat(paste("Observed data: ", DATA, "\n", sep=""))
cat(paste("Grid name: ", GRID, "\n", sep=""))
cat(paste("IDW power: ", IDWPOWER, "\n", sep=""))
cat(paste("Output file: ", OUTPUT, "\n", sep=""))

cat("\n")

############################################## Load packages

cat("Loading packages...\n")

load_pckg <- function(pkg = "raster"){
  suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
  cat(paste(pkg, packageVersion(pkg), "\n"))  
}

load_pckg("data.table")
load_pckg("sp")
load_pckg("sf")
load_pckg("raster")
load_pckg("gstat")

cat("\n")



############################################## Main workflow

## Main function
calc_resids <- function(datt, grid, idwpower = 1.4){
  # e.g.
  # datt <- "Residuals_ECM.RData"
  # grid <- "PredicitonChunks/Chunk_001.RData"
  # idwpower = 1.4

  ## Load points for residuals interpolation (one of the chunks)
  cat("Loading new prediction grid...\n")
  grid <- readRDS(grid)
  clz <- c("cellid", "LON", "LAT")
  grid <- grid[, ..clz]

  ## Convert to `Spatial` class
  gridd <- SpatialPointsDataFrame(
    coords = grid[, c("LON", "LAT")],
    data = grid,
    proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

  ## Load data residuals
  cat("Loading residuals...\n")
  datt <- readRDS(datt)

  ## Remove missing values
  if(any(is.na(datt$Residuals))){
    cat("Missing value removal...\n")
    datt <- datt[!is.na(datt$Residuals), ]
  }

  ## Convert to `Spatial` class
  xy <- SpatialPointsDataFrame(
    coords = datt[, c("LON", "LAT")],
    data = datt,
    proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

  ## Inverse distance weighted interpolation
  cat("Interpolating residuals to a grid...\n")
  KRG <- gstat::idw(
        formula = Residuals ~ 1,
        locations = xy,
        newdata = gridd,
        idp = idwpower
        )

  cat("Preparing results...\n")
  res <- data.table(
    coordinates(gridd),
    cellid = gridd$cellid,
    Residuals = KRG$var1.pred)

  return(res)
}


RES <- calc_resids(
  datt = DATA,
  grid = GRID,
  idwpower = IDWPOWER)

cat("\nSaving results...\n")
saveRDS(
  object = RES,
  file = OUTPUT,
  compress = "xz")


## Check time
end_time <- Sys.time()

tmm <- as.numeric(difftime(end_time, start_time, units = "min"))
cat("\nElapsed time: ", tmm, " minutes\n")


cat("\n")
cat("Session info:\n")
sessionInfo()
cat("\n")

