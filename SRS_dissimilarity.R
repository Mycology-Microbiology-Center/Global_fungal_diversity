#!/usr/bin/Rscript

## Usage:
# ./SRS_dissimilarity.R \
#    --input "ALLFUNGITABLE_phyloseq.RData" \
#    --depth 1000


############################################## Parse input parameters

## Check time
start_time <- Sys.time()


cat("Parsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse)) # don't say "Loading required package: optparse"

## Parse arguments
option_list <- list(
  make_option(c("-i", "--input"), action="store", default=NA, type='character', help="Input phyloseq"),
  make_option(c("-d", "--depth"), action="store", default=100L, type='integer', help="SRS depth")
  )
opt <- parse_args(OptionParser(option_list=option_list))


## Assign variables
INPUT <- opt$input
RARDEPTH <- as.integer( opt$depth )

OUTPUT <- paste0("SRS__", "Depth_", RARDEPTH)

## Log assigned variables
cat(paste("Input phyloseq: ", INPUT, "\n", sep=""))
cat(paste("SRS depth: ", RARDEPTH, "\n", sep=""))
cat(paste("Output file prefix: ", OUTPUT, "\n", sep=""))
cat("\n")

############################################## Load packages and data

cat("Loading packages...\n")

load_pckg <- function(pkg = "raster"){
  suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
  cat(paste(pkg, packageVersion(pkg), "\n"))
}

load_pckg("vegan")
load_pckg("phyloseq")
load_pckg("metagMisc")
load_pckg("SRS")
load_pckg("betapart")
load_pckg("plyr")

cat("\n")

set.seed(14789)

## Load snow and ice mask
cat("\nLoading input data..\n")
PHYS <- readRDS(INPUT)

############################################## Main workflow

## SRS-normalization (Scaling with Ranked Subsampling)
cat("SRS..\n")
RAR <- phyloseq_SRS(physeq = PHYS, Cmin = RARDEPTH,  drop_zeros = TRUE)

cat("Extracting OTU table..\n")
## betapart input = data frame, where rows are sites and columns are species
OTUTAB <- as.data.frame(otu_table(RAR))
if(taxa_are_rows(RAR) == TRUE){
  OTUTAB <- as.data.frame( t(OTUTAB) )
}

cat("Converting to presence-absence..\n")
PA <- decostand(x = OTUTAB, method = "pa")


cat("Estimating dissimilarities..\n")
DISM <- list()

## Incidence-based pair-wise dissimilarities
cat("..Sorenson\n")
DISM$sor <- try( beta.pair(PA, index.family = "sorensen") )

## Abundance-based pair-wise dissimilarities
cat("..Bray\n")
DISM$brc <- try( beta.pair.abund(OTUTAB, index.family = "bray") )

cat("Exporting results..\n")

saveRDS(object = RAR,
  file = paste0(OUTPUT, "_RAREFIED.RData"),
  compress = "xz")

saveRDS(object = DISM,
  file = paste0(OUTPUT, "_DISSIM.RData"),
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

