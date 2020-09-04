library(minfi, quietly=TRUE)
library(optparse, quietly=TRUE)
library(RColorBrewer, quietly=TRUE)
library(qs, quietly=TRUE)


# ---- 0. Parse CLI arguments
option_list <- list(
    make_option(c('-i', '--methylset'), 
        help=c("Path the the preprocessed MethylSet object "), 
        type="character"),
    make_option(c('-n', '--normalized'), 
        help=c("Path to the GenomicMethylSet containging the normalzied output from QC step 3."), 
        type='character'),
    make_option(c('-o', '--output'), 
        help='Path to write the pdf file containing the normalized vs unnormalized QC2 plots',
        type='character')
)

opt <- parse_args(OptionParser(option_list=option_list))


# ---- 1. Read in methylSet

methylSet <- qread(opt$methylset)


# ---- 2. Define QC Functions

diffMethylVsUnmethylPeaks <- function(betaValues) {
    # Drop any NA values
    betaValues <- na.omit(betaValues)

    # Calculate the density values
    densityResults <- density(betaValues)
    densityBetas <- densityList$x
    densityValues <- densityList$y

    # Find the middle index
    middleIndex <- floor(length(densityValues) / 2)

    # Find the index of the max density above and below the middle index
    whichMaxMethylated <- middleIndex + which.max(densityValues[seq(middleIndex + 1, length(densityValues))])
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

    # Range of beta values is (0, 1), therefore AUC of density plot is just number of obs over a specified interval
    methylAUC <- sum(betaValues >= 0.75)
    noiseAUC <- sum(betaValues < 0.75 & betaValues > 0.25)
    unmethylAUC <- sum(betaValues <= 0.25)

    signalOverNoise <- (methylAUC + unmethylAUC) / noiseAUC
    return(signalOverNoise)
}

# ---- 3. Calculate QC Metrics
