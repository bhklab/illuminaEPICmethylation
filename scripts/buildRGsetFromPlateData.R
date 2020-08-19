library(minfi)
library(optparse)
library(BiocParallel)
library(qs)


# ---- 0. Parse CLI arguments

option_list <- list(
    make_option(c('-p', '--plates'), 
        help=c("A comma separated list of paths to the directory for each microarray plate."), 
        type="character"),
    make_option(c('-l', '--labels'), 
        help=c("A comma separated list of paths to the labels for each plate, in the same order as '--plates'."), 
        type="character"),
    make_option(c('-o', '--output'), 
        help='Path and filename to save the output to.', 
        type='character'),
    make_option(c('-n', '--nthread'), 
        help=c("The number of threads to parallelize over."), 
        type='integer', default=1)
)

opt <- parse_args(OptionParser(option_list=option_list))


# ---- 1. Read in data
message("Reading in plate data...\n")


plateDirs <- unlist(strsplit(opt$plates, ' '))
arrays <- lapply(plateDirs, FUN=read.metharray.sheet)


# ---- 2. Label array data
message("Reading in plate labels...\n")


## FIXME:: Why is this making a list of lists?
labelPaths <- unlist(strsplit(opt$labels, ' '))

labels <- lapply(labelPaths, read.csv)

message("Labelling plate data...\n")

.labelArray <- function(array, label) {
    array$Sample_Group <- label$Group
    return(array)
}

labelledArrays <- mapply(FUN=.labelArray,  # Function to apply
                         array=arrays, label=labels,  # Arguments to map to (i.e., index N of all arguments passed to FUN)
                         SIMPLIFY=FALSE)  # Don't try to simplify to a matrix


# ---- 3. Read data into RGChannelSet objects
message("Building RGChannelSets from array data...\n")

rgSets <- bplapply(labelledArrays, function(targets) read.metharray.exp(targets=targets))


# --- 4. Annotate the RGChannelSet objects
message("Annotating RGChannelSets...\n")

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
message("Merging RGSets into one object...\n")

.combineArrays <- function(x1, x2) combineArrays(x1, x2, outType='IlluminaHumanMethylationEPIC')

finalRGSet <- Reduce(f=.combineArrays, finalRGSets)

# ---- 6. Save the RGSet to disk
message("Saving merged RGSet to disk...\n")

qsave(finalRGSet, file=opt$output)

message("Done!\n\n")