#!/usr/bin/env Rscript

## Dissimililarity summary

## Check time
start_time <- Sys.time()

cat("Parsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse)) # don't say "Loading required package: optparse"

## Parse arguments
option_list <- list(
  make_option(c("-d", "--dist"),   action="store", default=NA, type='character', help="Dissimilarity matrix"),
  make_option(c("-e", "--env"),    action="store", default=NA, type='character', help="Environmental data"),
  make_option(c("-o", "--output"), action="store", default=NA, type='character', help="Output suffix"),
  make_option(c("-g", "--geotransform"), action="store", default=TRUE, type='logical', help="Sqrt-transform geo distances")
  )
opt <- parse_args(OptionParser(option_list=option_list))


## Assign variables
dd  <- opt$dist
env <- opt$env
TRANSF <- opt$geotransform
OUTPUTSUFFIX <- opt$output

## Log assigned variables
cat(paste("Dissimilarity matrix: ", dd,    "\n", sep=""))
cat(paste("Environmental data: ",  env,    "\n", sep=""))
cat(paste("Geo-transformation: ",  TRANSF, "\n", sep=""))
cat(paste("Output suffix: ", OUTPUTSUFFIX, "\n", sep=""))

cat("\n")

############################################## Load packages and data

cat("Loading packages...\n")

load_pckg <- function(pkg = "raster"){
  suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
  cat(paste(pkg, packageVersion(pkg), "\n"))  
}

load_pckg("gdm")

cat("\n")


## Load pre-computed community dissimilarities
cat("Loading dissimilarities...\n")
dd <- readRDS(dd)

## Load environmental predictors
cat("Loading environmental predictors...\n")
env <- readRDS(env)

cat("Creating site-pair table\n")

## Convert pre-calculated `dist` to matix
dm <- as.data.frame( as.matrix(dd) )

## Check if all samples in dist are in env
if(any(!colnames(dm) %in% env$site)){
  warning("WARNING: some samples are missing from the env table\n")

  smps <- rownames(dm)
  smpsrm <- which(! smps %in% env$site)
  cat(".. Missing samples: ", paste(smps[ smpsrm ], collapse = " "), "\n")

  dm <- dm[-smpsrm, -smpsrm]
}

## Check if all samples in env are in dist
if(any(!env$site %in% colnames(dm))){
  warning("WARNING: some samples are missing from the dissim\n")
}

## Add a site ID column
site <- unique(colnames(dm))
dm <- cbind(site, dm)

## Reorder sample metadata
env <- env[ match(x = dm$site, table = env$site), ]

## Create GDM-formatted site-pair table
sitepair <- formatsitepair(
  bioData = dm,
  bioFormat = 3,
  predData = env,
  siteColumn = "site", XColumn = "LONG", YColumn = "LAT",
  verbose = TRUE)


cat(".. estimating Haversine distance\n")

geoDist <- geodist::geodist(
  x = env[, c("LONG", "LAT")],
  paired = FALSE, measure = "haversine")

geoDist <- data.frame(
  site = env$site,
  as.data.frame( geoDist / 1000 ))   # + convert dist to km

colnames(geoDist)[-1] <- env$site

cat(".. adding custom geographic distance to an existing site-pair table\n")

sitepair <- formatsitepair(
  bioData = sitepair,  # existing site-pair table
  bioFormat = 4,
  predData = env,      # should be provided, but will not be used
  siteColumn = "site",
  distPreds = list(geoDist)
  )

## Rename geodist columns
colnames(sitepair)[ colnames(sitepair) %in% "s1.matrix_1" ] <- "s1.GeoDist"
colnames(sitepair)[ colnames(sitepair) %in% "s2.matrix_1" ] <- "s2.GeoDist"


## Square-root transform geographic distance
## (Legendre et al., 2015, DOI:10.1111/2041-210X.12425)
if(TRANSF == TRANSF){
  cat(".. sqrt-transforming geographic distances\n")
  sitepair$s2.GeoDist <- sqrt(sitepair$s2.GeoDist)
}

cat(".. site-pair table created\n")

## Remove samples with missing data
if(any(is.na(sitepair))){
  cat("..Removing samples with missing data\n")
  sitepair <- na.omit(sitepair)
}

## Export site-pair table
cat("..Exporting site-pair table\n")
saveRDS(object = sitepair,
  file = paste0("GDM_SitePair_", OUTPUTSUFFIX, ".RData"),
  compress = "xz")


cat("Fitting a Generalized Dissimilarity Model\n")

gdmMod <- gdm(sitepair,
  geo = FALSE,     # geo is taken from the custom distance matrix (GeoDist)
  splines = NULL,  # splines = rep(3, times = length(colnames(sitepair)[-c(1:6)])/2)
  knots = NULL)

cat("..done\n")

# summary(gdmMod)
# plot(gdmMod, plot.layout = c(1,2))

cat("Exporting model\n")

saveRDS(object = gdmMod,
  file = paste0("GDM_model_", OUTPUTSUFFIX, ".RData"),
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

