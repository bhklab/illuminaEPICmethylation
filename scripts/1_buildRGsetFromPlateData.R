# ---- 0. Load dependencies
message("Loading script dependencies...\n")

# Suppress package load messages to ensure logs are not cluttered
suppressMessages({
    library(minfi, quietly=TRUE)
    library(optparse, quietly=TRUE)
    library(BiocParallel, quietly=TRUE)
    library(qs, quietly=TRUE)
})

# ---- 0. Parse CLI arguments

input <- snakemake@input
output <- snakemake@output

cat(input)
cat(output)

# ---- 1. Read in data
message("Reading in plate data...\n")


plateDirs <- unlist(strsplit(input$plates, ' '))
arrays <- lapply(plateDirs, FUN=read.metharray.sheet)


# ---- 2. Label array data
message("Reading in plate labels...\n")


## FIXME:: Why is this making a list of lists?
labelPaths <- unlist(strsplit(input$labels, ' '))

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

rgSets <- bplapply(labelledArrays, function(targets) read.metharray.exp(targets=targets, extended=TRUE))


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

.combineArrays <- function(rgSet1, rgSet2) {
    combineArrays(rgSet1, rgSet2, outType='IlluminaHumanMethylationEPIC')
}

finalRGSet <- Reduce(f=.combineArrays, finalRGSets)

# ---- 6. Save the RGSet to disk
message("Saving merged RGSet to disk...\n")

qsave(finalRGSet, file=output[[1]])

message("Done!\n\n")