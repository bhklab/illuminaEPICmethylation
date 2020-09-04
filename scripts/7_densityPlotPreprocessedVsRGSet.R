# ---- 0. Load dependencies
message('Loading script dependencies...\n')

# Suppress package load messages to ensure logs are not cluttered
suppressMessages({
    library(minfi, quietly=TRUE)
    library(qs, quietly=TRUE)
    library(optparse, quietly=TRUE)
})

# ---- 0. Parse CLI arguments
message('Parsing command line arguments...\n')

option_list <- list(
    make_option(c('-m', '--methylsets'), 
        help=c("Path the the preprocessed MethylSet object "), 
        type="character"),
    make_option(c('-r', '--rgset'), 
        help=c("Path to the RGSets from the previous two QC steps"), 
        type='character'),
    make_option(c('-o', '--output'), 
        help='Path to write the pdf of density plots for each MethylSet.',
        type='character'),
    make_option(c('-t', '--titles'),
        help='Title on plots for each methylSet included in -m.',
        type='character')
)

opt <- parse_args(OptionParser(option_list=option_list))

methylSetPaths <- unlist(strsplit(opt$methylsets, split=' '))
plotTitles <- c('Unprocessed', unlist(strsplit(opt$titles, split=' ')))

# ---- 1. Read in methylSet
message("Reading in data from:\n\t", 
        paste0(opt$rgset, '\n\t', paste0(methylSetPaths, collapse='\n\t'), '\n'))

methylSets <- lapply(methylSetPaths, qread)
rgSet <- qread(opt$rgset)

plotSets <- c(list(rgSet), methylSets)

# ---- 2. Define QC Functions

diffMethylVsUnmethylPeaks <- function(betaValues) {
    # Drop any NA values
    betaValues <- na.omit(betaValues)

    # Calculate the density values
    densityResults <- density(betaValues)
    densityBetas <- densityResults$x
    densityValues <- densityResults$y

    # Find the middle index
    middleIndex <- floor(length(densityValues) / 2)

    # Find the index of the max density above and below the middle index
    whichMaxMethylated <- middleIndex + 
        which.max(densityValues[seq(middleIndex + 1, length(densityValues))])
    whichMaxUnmethylated <- which.max(densityValues[seq_len(middleIndex)])

    # Get the peak of methylated vs unmethylated signal
    methylatedPeak <- densityBetas[whichMaxMethylated]
    unmethylatedPeak <- densityBetas[whichMaxUnmethylated]

    peakDifference <- methylatedPeak - unmethylatedPeak
    return(peakDifference)
}

calcMethylSignalOverNoise <- function(betaValues) {
    # Drop any NA values
    betaValues <- na.omit(betaValues)
    
    # Range of beta values is (0, 1), therefore AUC of density plot is 
    #   just number of obs over a specified interval
    methylAUC <- sum(betaValues >= 0.75)
    noiseAUC <- sum(betaValues < 0.75 & betaValues > 0.25)
    unmethylAUC <- sum(betaValues <= 0.25)

    signalOverNoise <- (methylAUC + unmethylAUC) / noiseAUC
    return(signalOverNoise)
}

# ---- 3. Calculate QC Metrics and Density Plot All Preprocess Methods for Each Sample

# Write functions to apply over each betaValue matrix
.sampleMethylVsUnmethylPeaks <- function(betaMatrix) {
    apply(betaMatrix, 2, FUN=diffMethylVsUnmethylPeaks)
}

.sampleMethylSignalOverNosie <- function(betaMatrix) {
    apply(betaMatrix, 2, FUN=calcMethylSignalOverNoise)
}

message("Extracting beta value matrices...\n")
# NOTE: This assumes the same samples for all plotSets
sampleNames <- colnames(plotSets[[1]])


.getBetaNaOmit <- function(methylSet) na.omit(getBeta(methylSet))
betaMatrixList <- lapply(plotSets, FUN=.getBetaNaOmit)

message("Calculating QC metrics...\n")
peakDifferences <- lapply(betaMatrixList, FUN=.sampleMethylVsUnmethylPeaks)
aucRatios <- lapply(betaMatrixList, FUN=.sampleMethylSignalOverNosie)

numRows <- ceiling(length(plotSets) / 2)

message(paste0("Writing plots to:\n\t", opt$output))
pdf(opt$output)

for (sampleIdx in seq_along(sampleNames)) {
    
    message(paste0('Plotting: ', sampleNames[sampleIdx]))
    print(sampleNames[sampleIdx])
    plot.new()
    par(mfrow=c(numRows, 2))
    
    for (plotIdx in seq_along(plotSets)) {
        plot(density(betaMatrixList[[plotIdx]][, sampleIdx]),
             title=plotTitles[plotIdx])
        mtext(paste0('Peak Diff: ', peakDifferences[[plotIdx]][sampleIdx]), 
              side=4, line=0)
        mtext(paste0('AUC Ratio: ', aucRatios[[plotIdx]][sampleIdx]),
              side=4, line=1)
    }
}

dev.off()

message("Done!\n\n")