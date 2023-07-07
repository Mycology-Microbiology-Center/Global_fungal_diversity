#!/usr/bin/env Rscript

## GDM - Estimate local turnover
## average compositional similarity within a NNN-m radius

## Example:
# ./GDM_LocalTurnover_script.R \
#   --raster    "GDMsm_SimpsNoAbund.tif" \
#   --spatvect  "SpatVectors/SpatVect_00001.gpkg" \
#   --buffer 150000 \
#   --threads   1


## Check time
start_time <- Sys.time()

cat("Parsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse)) # don't say "Loading required package: optparse"

## Parse arguments
option_list <- list(
  make_option(c("-r", "--raster"),   action="store", default=NA,   type='character', help="GDM raster (non-scaled PCA scores)"),
  make_option(c("-s", "--spatvect"), action="store", default=NA,   type='character', help="Point coordinates (GeoPackage format)"),
  make_option(c("-b", "--buffer"),   action="store", default=150000, type='integer', help="Spatial buffer size (in meters)"),
  make_option(c("-t", "--threads"),  action="store", default=8L, type='integer', help="Number of CPUs to use"),
  make_option(c("-o", "--output"),   action="store", default=NA, type='character', help="Output suffix")
  )
opt <- parse_args(OptionParser(option_list=option_list))

## Assign variables
RASTER   <- opt$raster
SPATVECT <- opt$spatvect
BUFFER   <- as.integer( opt$buffer )
THREADS  <- as.integer( opt$threads )
OUTPUTSUFFIX <- opt$output

OUTPUTNAME <- paste0(OUTPUTSUFFIX, "_",
  sub(pattern = ".gpkg", replacement = "", x = basename(SPATVECT)))

## Log assigned variables
cat("\n")
cat(paste("GDM input raster: ",   RASTER,   "\n", sep=""))
cat(paste("Points of interest: ", SPATVECT, "\n", sep=""))
cat(paste("Buffer size: ",        BUFFER,   "\n", sep=""))
cat(paste("Number of CPU threads: ", THREADS, "\n", sep=""))
cat(paste("Output suffix: ",    OUTPUTSUFFIX, "\n", sep=""))
cat(paste("Output file name: ", OUTPUTNAME,   "\n", sep=""))

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

if(THREADS > 1){
  cat("Setting up local cluster\n")

  load_pckg("doFuture")
  registerDoFuture()
  # plan(multisession, workers = 10)     # for RStudio
  plan(multicore, workers = THREADS)     # will crash RStudio
  options(future.globals.maxSize = 5e9)  # 5GB; default = 500 * 1024 ^ 2 = 500 MiB

}



########################################## Main workflow

## Load raster
cat("Loading raster...\n")
RST <- rast(RASTER)

## Load data
cat("Loading spatial points...\n")
dtt <- vect(SPATVECT)

## Main function
local_dissim <- function(x, bufsize_meters = BUFFER){
  # x = single point of class `SpatVector` (terra); e.g., x <- dtt[1,]
  # bufsize_meters = buffer size in meters, default = 150 km; e.g. bufsize_meters <- 150000

  ## Add buffer around the point
  pl <- terra::buffer(x = x, width = bufsize_meters)

  ## Extract values from SpatRaster by SpatVector
  rr <- extract(x = RST, y = pl, na.rm = TRUE)
  rr <- na.omit(rr)

  ## Estimate and summarize distances
  dd <- dist(rr[, 2:4])
  res <- data.table(
    cellid = x$cellid,
    DistMedian = median(dd, na.rm = TRUE),
    DistQ1 = quantile(dd, probs = 0.25, na.rm = TRUE),
    DistQ3 = quantile(dd, probs = 0.75, na.rm = TRUE)
    )
  rm(dd, pl, rr)
  return(res)
}

cat("\nRunning analysis\n")

## Process all points
RES <- alply(
  .data = 1:nrow(dtt),
  .margins = 1,
  .fun = function(x){
    res <- try( local_dissim(dtt[x, ]) )
    if("try-error" %in% class(res)){
      res <- data.table(cellid = dtt[x, ]$cellid,
        DistMedian = NA, DistQ1 = NA, DistQ3 = NA)
    }
    return(res)
  }, .parallel = TRUE)

RES <- rbindlist(RES)


cat("Exporting results\n")

saveRDS(object = RES,
  file = paste0(OUTPUTNAME, ".RData"),
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
