library(minfi)
library(optparse)
library(qs)


# ---- 0. Parse CLI arguments

option_list <- list(
    make_option(c('-p', '--plates'), 
        help=c("A comma separated list of paths to the directory for each microarray plate."), 
        type="character"),
    make_option(c('-l', '--labels'), 
        help=c("A comma separated list of paths to the labels for each plate, in the same order as '--plates'."), 
        type="character"),
    make_option(c('-n', '--nthread'), 
        help=c("The number of threads to parallelize over.", 
        type='integer', default=1))
)

opts <- OptionParser(option_list=option_list)


# ---- 1. Read in data
message("Reading in plate data...")

plateDirs <- strsplit(opts$plates, ',')
arrays <- lapply(plateDirs, FUN=read.metharray.sheet)


# ---- 2. Label array data
message("Reading in plate labels...")

labelsPaths <- strsplit(opts$labels, ',')
labels <- lapply(labelPaths, FUN=read.csv, stringsAsFactors=FALSE)

message("Labelling plate data...")

.labelArray <- function(array, label) {
    array$Sample_Group <- label$Group
    return(array)
}

labelledArrays <- mapply(FUN=.labelArray,  # Function to apply
                         array=arrays, label=labels,  # Arguments to map to (i.e., index N of all arguments passed to FUN)
                         SIMPLIFY=FALSE)  # Don't try to simplify to a matrix


# ---- 3. Read data into RGChannelSet objects
message("Building RGChannelSets from array data...")

rgSets <- lapply(labelledArrays, FUN=read.metharray.exp)


# --- 4. Annotate the RGChannelSet objects
message("Annotating RGChannelSets...")

.buildSampleName <- function(labelledArray) {
    paste(labelledArray$Sample_Plate, labelledArray$Sample_Well, labelledArray$Sample_Group, sep='-')
}

sampleLabels <- lapply(labelledArrays, FUN=.buildSampleName)

.labelRGSet <- function(rgSet, sampleLabel) { 
    sampleNames(rgSet) <- sampleLabel
    return(rgSet)
}

labelledRGSets <- mapply(FUN=.labelRGSet, 
                         rgSet=rgSets, sampleLabel=sampleLabels,
                         SIMPLIFY=FALSE)

.addMetadata <- function(rgSet, metadata) {
    rgSet$metadata  <- metadata
    return(rgSet)
}

finalRGSets <- mapply(FUN=.addMetadata,
                      rgSet=labelledRGSets, metadata=labels,
                      SIMPLIFY=FALSE)


# ---- 5. Merge the RGSets into a single object
