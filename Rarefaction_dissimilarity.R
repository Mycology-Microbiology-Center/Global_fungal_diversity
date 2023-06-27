#!/usr/bin/Rscript
##  #!/gpfs/space/home/kupagme/miniconda3/bin/Rscript

## Usage:
# ./Rarefaction_dissimilarity.R \
#    --input "ALLFUNGITABLE_phyloseq.RData" \
#    --depth 1000 \
#    --seed 1



############################################## Parse input parameters

## Check time
start_time <- Sys.time()


cat("Parsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse)) # don't say "Loading required package: optparse"

## Parse arguments
option_list <- list(
  make_option(c("-i", "--input"), action="store", default=NA, type='character', help="Input phyloseq"),
  make_option(c("-d", "--depth"), action="store", default=100L, type='integer', help="Rarefaction depth"),
  make_option(c("-s", "--seed"), action="store", default=1L, type='integer', help="Random seed")
)
opt <- parse_args(OptionParser(option_list=option_list))


## Assign variables
INPUT <- opt$input
RARDEPTH <- as.integer( opt$depth )
SEED <- as.integer( opt$seed )

OUTPUT <- paste0("Depth_", RARDEPTH, "_Seed_", SEED)

## Log assigned variables
cat(paste("Input phyloseq: ", INPUT, "\n", sep=""))
cat(paste("Rarefaction depth: ", RARDEPTH, "\n", sep=""))
cat(paste("Random seed: ", SEED, "\n", sep=""))
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
# load_pckg("metagMisc")
load_pckg("betapart")
load_pckg("plyr")

cat("\n")

set.seed(SEED)

## Load snow and ice mask
cat("\nLoading input data..\n")
PHYS <- readRDS(INPUT)

############################################## Main workflow

cat("Rarefying..\n")
RAR <- phyloseq::rarefy_even_depth(PHYS,
  sample.size = RARDEPTH, rngseed = SEED, 
  trimOTUs = TRUE, replace = FALSE,
  verbose = TRUE)


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
# cat("..Jaccard\n")
# DISM$jac <- try( beta.pair(PA, index.family = "jaccard") )

## Abundance-based pair-wise dissimilarities
cat("..Bray\n")
DISM$brc <- try( beta.pair.abund(OTUTAB, index.family = "bray") )
# cat("..Ruzicka\n")
# DISM$rzk <- try( beta.pair.abund(OTUTAB, index.family = "ruzicka") )


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

