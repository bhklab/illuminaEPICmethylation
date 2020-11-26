# ---- 0. Load dependencies
message('Loading script dependencies...\n')

# Suppress package load messages to ensure logs are not cluttered
suppressMessages({
    library(minfi, quietly=TRUE)
    library(data.table, quietly=TRUE)
    library(qs, quietly=TRUE)
    library(wateRmelon, quietly=TRUE)
})

# ---- 0. Parse Snakemake arguments
input <- snakemake@input
params <- snakemake@params
output <- snakemake@output


# ---- 1. Read in rgSet
message(paste0('Reading input file from ', input$rgset, '...\n'))

rgSet <- qread(input$rgset)


methods <- params$preprocess_methods
outputs <- output$methylsets
reports <- output$qc_reports

for (i in seq_along(outputs)) {

    # ---- 2. Preprocess the rgSet to one or more MethylSets using the selected preprocessing methods
    method <- paste0('preprocess', methods[i])
    message(paste0("Preprocessing with ", method, " function...\n"))

    .preprocessFunction <- get(method)  # Gets the object named the same as the string passed to it

    if (method == 'preprocessFunnorm') {
        methylSet <- .preprocessFunction(rgSet, ratioConvert=FALSE)
    } else {
        methylSet <- .preprocessFunction(rgSet)
    }

    # ---- 3. Run QC on the methylSet
    message('Running minfiQC on the methylSet...\n')

    qcResults <- minfiQC(methylSet, verbose=TRUE)


    # ---- 4. Extract the resulting methylSet and QC file
    message("Extracting QC results...\n")

    methylSet <- qcResults[[1]]

    qcDf <- qcResults[[2]]


    # ---- 5. Write output files
    message(paste0("Saving methylSet to ", outputs[i], '...\n'))

    qsave(methylSet, file=outputs[i])

    message(paste0("Saving QC results to .csv in ", reports[i], '...\n'))
    write.csv(qcDf, file=reports[i])

    message("\n\n")
}


message("Done!\n\n")