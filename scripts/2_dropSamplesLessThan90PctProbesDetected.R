# ---- 0. Load dependencies
message("Loading script dependencies...\n")

# Suppress package load messages to ensure logs are not cluttered
suppressMessages({
    library(minfi, quietly=TRUE)
    library(IlluminaHumanMethylationEPICmanifest, quietly=TRUE)
    library(data.table, quietly=TRUE)
    library(qs, quietly=TRUE)
})

# ---- 0. Parse Snakemake arguments

input <- snakemake@input
params <- snakemake@params
output <- snakemake@output


# ---- 1. Load RGChannelSet
message("Reading in RGChannelSet from: ", input$rgset, '...\n')
rgSet <- qread(input$rgset)


# ---- 2. Calculate detection p-value, proportion/count failed probes per sample, probes with proportion failed samples 
message("Calculating QC metrics...\n")
source('scripts/functions/summarizeDetectionPvalueQC.R')

detectionQCResults <- summarizeDetectionPvalueQC(rgSet, pValue=params$detection_pvalue)


# ---- 4. Construct result names and save QC results to disk

detectionPOut <- output$detectionPValues
sampleQCOut <- output$sampleQC
probeQCOut <- output$probeQC

paths <- c(detectionPOut, sampleQCOut, probeQCOut)

for (i in seq_along(detectionQCResults)) {
    message(paste0('\tSaving QC results to: ', paths[i], '...\n'))
    fwrite(detectionQCResults[[i]], file=paths[i])
}


# ---- 5. Filter failed samples from RGChannelSet and save to disk
message(paste0('Filtering RGChannelSet to samples with >90% CpGs detected at alpha ', params$detection_pvalue, '...\n'))
# Drop rownames column so all values are numeric and we can use colMeans
detPvals <- detectionQCResults$detectionPvals[, -'rownames']
detected <- detPvals < params$detection_pvalue


numSamplesPassedQC <- sum(colMeans(detected) > 0.9)
message(paste0(numSamplesPassedQC, 'out of ', ncol(detPvals), ' samples passed QC at alpha of ', params$detection_pvalue, '...\n'))

failedSamples <- which(colMeans(detected) < 0.9)
message(paste0("The following samples failed QC: \n\t", paste0(names(failedSamples), collapse=',\n\t'), '\n'))

keepSamples <- colMeans(detected) > 0.9

rgSetFiltered <- rgSet[, keepSamples]
print(rgSetFiltered)

message(paste0('Saving filtered RGSet to: ', output$rgset_filtered, '...\n'))
qsave(rgSetFiltered, file=output$rgset_filtered)


message("Done!\n\n")