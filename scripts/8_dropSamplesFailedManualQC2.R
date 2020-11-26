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

input <- snakemake@input
output <- snakemake@output
params <- snakemake@params

# ---- 1. Read in Preprocessed MethylSet
message(paste0("Loading MethylSet object from ", input$methylset, '...\n'))
methylSet <- qread(input$methylset)

message(paste0("Loading QC Metrics from ", input$qc_report, '...\n'))
qcMetricDT <- fread(input$qc_report)[preprocMethod == params$preproc_method, ]

# ---- 2. Subset out failed samples
message('Removing failed samples...\n')

for (i in seq_along(params$qc_criteria)) {

    print(names(params$qc_criteria)[i])
    switch(names(params$qc_criteria)[i],
        'minAucRatio'={ keepSamples <- qcMetricDT[aucRatio > params$qc_criteria$minAucRatio, ]$sample },
        'maxNumPeaks'={ keepSamples <- qcMetricDT[numPeaks < params$qc_criteria$maxNumPeaks, ]$sample },
        'minPeakDiff'={ keepSamples <- qcMetricDT[peakDiff > params$qc_criteria$minPeakDiff, ]$sample },
        stop("Unknown QC filtering criteria, valid options are minAucRatio, maxNumPeaks 
            an maxPeakDiff"))

    print(keepSamples)

    # Note: need to make.names because we formatted the sample names for plotting in the previous rule
    methylSet <- methylSet[, make.names(colnames(methylSet)) %in% keepSamples]
    qcMetricDT <- qcMetricDT[sample %in% keepSamples]

    qsave(methylSet, file=output[[i]])
}

message("Done!\n\n")
