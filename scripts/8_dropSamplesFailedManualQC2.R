# ---- 0. Load dependencies
message("Loading script dependencies...\n")

# Suppress package load messages to ensure logs are not cluttered
suppressMessages({
    library(minfi, quietly=TRUE)
    library(data.table, quietly=TRUE)
    library(qs, quietly=TRUE)
    library(optparse, quietly=TRUE)
})

# ---- 0. Parse CLI arguments
message("Parsing command line arguments...\n")
option_list <- list(
    make_option(c('-m', '--methylset'), 
        help=c("Path to the input RGChannelSet object."), 
        type="character"),
    make_option(c('-q', '--qc_steps'),
        help='String specifying the manual QC steps to apply, 
              as specified in the manual_qc_steps field of config.yml',
        type='character'),
    make_option(c('-o', '--output'),
        help='Path and filename to save the per sample probe QC results to.', 
        type='character')
)

opt <- parse_args(OptionParser(option_list=option_list))

output <- unlist(strsplit(opt$output, split=' '))

# ---- 1. Read in RGSet
message(paste0("Loading RGChannelSet object from ", opt$rgset, '...\n'))
rgSet <- qread(opt$rgset)


# ---- 2. Subset out failed samples
message('Removing failed samples...\n')

colData <- data.table(as.data.frame(colData(rgSet)))

.strsplitUnlist <- function(str, split) unlist(strsplit(str, split=split), recursive=FALSE)

input <- .strsplitUnlist(opt$qc_steps, split=' ')
steps <- lapply(input, FUN=.strsplitUnlist, split=':')
manualQCStepNames <- unlist(lapply(steps, `[[`, i=1))
samplesFailedQCPerStep <- lapply(steps, `[[`, i=2)

samplesFailedQCPerStepPerPlate <- lapply(samplesFailedQCPerStep, .strsplitUnlist, split=';')

.strsplitList <- function(list) {
    lapply(list, .strsplitUnlist, split=',')
}

# NOTE: This assumes that the sorted plates are in the same order as samples/wells 
#   in the config.yml file!
failedSamplesList <- lapply(samplesFailedQCPerStepPerPlate, 
                            FUN=.strsplitList)

allPlateNames <- sort(unique(colData$Sample_Plate))

.getPlateWellIndexes <- function(plate, wells, colData) {
    return(colData[Sample_Plate == plate & Sample_Well %in% wells, which=TRUE])
}

for (i in seq_along(failedSamplesList)) {
    message(paste0('Failed samples for manual qc step: ', manualQCStepNames[i], '...\n'))
    print(failedSamplesList[[i]])
    plateNames <- allPlateNames[seq_along(failedSamplesList[[i]])]

    dropIndexes <- unlist(mapply(
        FUN=.getPlateWellIndexes,  # Function to apply
        plate=plateNames, wells=failedSamplesList[[i]],  # Arguments to map over together
        MoreArgs=list('colData'=colData),  # Additional arguments
        SIMPLIFY=FALSE  # Don't try to simplify to a matrix or vector
        ))

    keepIndexes <- setdiff(seq_len(ncol(rgSet)), dropIndexes)

    message("Removing samples failed QC...\n")
    rgSet <- rgSet[, keepIndexes]
    print(rgSet)
    message("\n")

    message(paste0('Saving RGChannelSet for ', manualQCStepNames[i], ' to: ', output[i], '...\n'))
    qsave(rgSet, output[i])
}


message("Done!\n\n")