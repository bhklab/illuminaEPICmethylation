library(minfi)
library(qs)
library(optparse)
library(matrixStats)
library(IlluminaHumanMethylationEPICmanifest)
library(data.table)


# ---- 0. Parse CLI arguments
option_list <- list(
    make_option(c('-i', '--input'), 
        help=c("Path to the input RGChannelSet object."), 
        type="character"),
    make_option(c('-d', '--detection'),
        help='Path and filename to detection p-value dataframe to.', 
        type='character'),
    make_option(c('-p', '--probes'),
        help='Path and filename to save the per sample probe QC results to.', 
        type='character'),
    make_option(c('-s', '--samples'),
        help='Path and filename to save the number of samples failed at various p-value cut-offs', 
        type='character'),
    make_option(c('-r', '--report'),
        help='Path and filename to save the minfi QC report .pdf to.', 
        type='character')
)

opt <- parse_args(OptionParser(option_list=option_list))


# ---- 1. Read in the RGSet from the previous step
message(paste0("Loading RGChannelSet object from", opt$input, '...\n'))

rgSet <- qread(opt$input)


# ---- 2. Summarize detection p-value qc metrics
message("Summarizing QC metrics based on detection p-value")

source("scripts/functions/summarizeDetectionPvalueQC.R")

detectionPvalueSummary <- summarizeDetectionPvalueQC(rgSet, pValue=0.01)

print(detectionPvalueSummary$sampleQC)
print(t(detectionPvalueSummary$probeQC))

# ---- 3. Write detection p-value qc metrics to .csv
message(paste0("Writing detection p-value qc metrics to: ", 
        paste0(opt$detection, opt$probes, opt$samples, collapse=', ')))

fwrite(detectionPvalueSummary$detectionPvals, file=opt$detection)
fwrite(detectionPvalueSummary$sampleQC, file=opt$samples)
fwrite(detectionPvalueSummary$probeQC, file=opt$probes)


# ---- 5. Generate a QC report for the RGset
message(paste0("Writing minfi QC report to ", opt$pdf, '...\n'))

qcReport(rgSet,
    sampNames=sampleNames(rgSet),
    sampGroups=colData(rgSet)$Sample_Group,
    pdf=opt$report
)


message("Done!\n\n")