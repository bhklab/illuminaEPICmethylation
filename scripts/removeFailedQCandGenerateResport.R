library(minfi)
library(IlluminaHumanMethylationEPICmanifest)
library(data.table)
library(qs)
library(optparse)


# ---- 0. Parse CLI arguments
option_list <- list(
    make_option(c('-i', '--input'), 
        help=c("Path to the input RGChannelSet object."), 
        type="character"),
    make_option(c('-f', '--failed'),
        help='Per plate failed samples', 
        type='character'),
    make_option(c('-d', '--detection'),
        help='Path and filename to detection p-value dataframe to.', 
        type='character'),
    make_option(c('-p', '--probes'),
        help='Path and filename to save the per sample probe QC results to.', 
        type='character'),
    make_option(c('-s', '--samples'),
        help='Path and filename to save the number of samples failed at various p-value cut-offs', 
        type='character'),
    make_option(c('-o', '--output'),
        help='Path and filename to save the per sample probe QC results to.', 
        type='character')
)

opt <- parse_args(OptionParser(option_list=option_list))


# ---- 1. Read in RGSet
message(paste0("Loading RGChannelSet object from ", opt$input, '...\n'))

rgSet <- qread(opt$input)


# ---- 2. Subset out failed samples
message('Removing failed samples...\n')

colData <- colData(rgSet)

perPlateFailedSamples <- strsplit(opt$failed, split=' ')
failedSampleList <- lapply(perPlateFailedSamples, FUN=strsplit, split=',')

# NOTE: This assumes that the sorted sample are in the same order as the config.yml file!
plateNames <- sort(unique(colData$Sample_Plates))
plateNames <- plateNames[seq_along(failedSampleList)]

.matchPlateAndWell <- function(colData, plateName, wellName) {
    colData$Sample_Plate == plateName & wellName %in% colData$Sample_Well 
}

dropSamples <- mapply(f=.matchPlateAndWell, plateName=plateNames, )


# ---- 3. Summarize detection p-value qc metrics
message("Summarizing QC metrics based on detection p-value...\n")

source('scripts/functions/summarizeDetectionPvalueQC.R')

detectionPvalueSummary <- summarizeDetectionPvalueQC(rgSet, pValue=0.01)

print(detectionPvalueSummary$sampleQC)
print(t(detectionPvalueSummary$probeQC))


# ---- 4. Write summarized detection p-value qc metrics to csv

message(paste0("Writing detection p-value qc metrics to: ", 
        paste0(opt$detection, opt$probes, opt$samples, collapse=', '), '\n'))

fwrite(detectionPvalueSummary$detectionPvals, file=opt$detection)
fwrite(detectionPvalueSummary$sampleQC, file=opt$samples)
fwrite(detectionPvalueSummary$probeQC, file=opt$probes)

message("Done!\n\n")