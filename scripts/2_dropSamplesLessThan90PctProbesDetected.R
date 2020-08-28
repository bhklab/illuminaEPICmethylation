# ---- 0. Load dependencies
message("Loading script dependencies...\n")

# Suppress package load messages to ensure logs are not cluttered
suppressMessages({
    library(minfi, quietly=TRUE)
    library(IlluminaHumanMethylationEPICmanifest, quietly=TRUE)
    library(data.table, quietly=TRUE)
    library(qs, quietly=TRUE)
    library(optparse, quietly=TRUE)
})


# ---- 0. Parse CLI arguments
option_list <- list(
    make_option(c('-r', '--rgset'), 
        help="Path to the input RGChannelSet object.", 
        type="character"),
    make_option(c('-p', '--pvalue'),
        help='P-value cut-off for calling probes as detected.',
        type='numeric'),
    make_option(c('-P', '--path'),
        help='Path to save the output files to.',
        type='character'),
    make_option(c('-o', '--output'),
        help='Path and filename to save the filtered RGChannelSet to.', 
        type='character')
)

opt <- parse_args(OptionParser(option_list=option_list))


# ---- 1. Load RGChannelSet
message("Reading in RGChannelSet from: ", opt$rgset, '...\n')
rgSet <- qread(opt$rgset)


# ---- 2. Calculate detection p-value, proportion/count failed probes per sample, probes with proportion failed samples 
message("Calculating QC metrics...\n")
source('scripts/functions/summarizeDetectionPvalueQC.R')

detectionQCResults <- summarizeDetectionPvalueQC(rgSet, pValue=opt$pvalue)


# ---- 4. Construct result names and save QC results to disk

detectionPOut <- paste0(opt$path, '.detection_pvalues.csv')
sampleQCOut <- paste0(opt$path, '.probes_failed_per_sample_p', opt$pvalue, '.csv')
probeQCOut <- paste0(opt$path, '.num_probes_with_proportion_failed_samples_p', opt$pvalue, '.csv')

paths <- c(detectionPOut, sampleQCOut, probeQCOut)

for (i in seq_along(detectionQCResults)) {
    message(paste0('\tSaving QC results to: ', paths[i], '...\n'))
    fwrite(detectionQCResults[[i]], file=paths[i])
}


# ---- 5. Filter failed samples from RGChannelSet and save to disk
message(paste0('Filtering RGChannelSet to samples with >90% CpGs detected at alpha ', opt$pvalue, '...\n'))
# Drop rownames column so all values are numeric and we can use colMeans
detPvals <- detectionQCResults$detectionPvals[, -'rownames']
detected <- detPvals < opt$pvalue


numSamplesPassedQC <- sum(colMeans(detected) > 0.9)
message(paste0(numSamplesPassedQC, 'out of ', ncol(detPvals), ' samples passed QC at alpha of ', opt$pvalue, '...\n'))

failedSamples <- which(colMeans(detected) < 0.9)
message(paste0("The following samples failed QC: \n\t", paste0(names(failedSamples), collapse=',\n\t'), '\n'))

keepSamples <- colMeans(detected) > 0.9

rgSetFiltered <- rgSet[, keepSamples]
print(rgSetFiltered)

message(paste0('Saving filtered RGSet to: ', opt$output, '...\n'))
qsave(rgSetFiltered, file=opt$output)


message("Done!\n\n")