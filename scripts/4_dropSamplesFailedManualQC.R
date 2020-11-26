# ---- 0. Load dependencies
message("Loading script dependencies...\n")

# Suppress package load messages to ensure logs are not cluttered
suppressMessages({
    library(minfi, quietly=TRUE)
    library(data.table, quietly=TRUE)
    library(qs, quietly=TRUE)
})

# ---- 0. Parse Snakemake arguments
input <- snakemake@input
params <- snakemake@params
output <- snakemake@output


print(params)

# ---- 1. Read in RGSet
message(paste0("Loading RGChannelSet object from ", input$rgset, '...\n'))
rgSet <- qread(input$rgset)


# ---- 2. Subset out failed samples
message('Removing failed samples...\n')

colData <- data.table(as.data.frame(colData(rgSet)))

.getPlateWellIndexes <- function(plate, wells, colData) {
    return(colData[Sample_Plate == plate & Sample_Well %in% wells, which=TRUE])
}

manual_qc_steps <- params$manual_qc_steps
qc_step_names <- names(manual_qc_steps)

for (i in seq_along(manual_qc_steps)) {
    message(paste0('Failed samples for manual qc step: ', qc_step_names[i], '...\n'))
    print(manual_qc_steps[[i]])
    plateNames <- names(manual_qc_steps[[i]])

    dropIndexes <- unlist(mapply(
        FUN=.getPlateWellIndexes,  # Function to apply
        plate=plateNames, wells=manual_qc_steps[[i]],  # Arguments to map over together
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
