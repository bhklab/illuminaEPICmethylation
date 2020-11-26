# ---- 0. Load dependencies
message('Loading script dependencies...\n')

# Suppress package load messages to ensure logs are not cluttered
suppressMessages({
    library(minfi, quietly=TRUE)
    library(data.table, quietly=TRUE)
    library(qs, quietly=TRUE)
    library(optparse, quietly=TRUE)
    library(wateRmelon, quietly=TRUE)
})


# ---- 0. Parse CLI arguments

option_list <- list(
    make_option(c('-i', '--rgset'),
        help=c("Path to the quality controlled RGChannelSet object."),
        type="character"),
    make_option(c('-m', '--methods'),
        help='Comma separated string specify the preprocessing methods to use for making MethylSets. 
            Options are: Illumina, SWAN, Noob and Funnorm.',
        type='character'),
    make_option(c('-o', '--outputs'),
        help=c("Path(s) to save the resulting methylSet(s) to."),
        type="character"),
    make_option(c('-r', '--reports'), 
        help='Path(s) to save the qc file(s) to.', 
        type='character')
)

opt <- parse_args(OptionParser(option_list=option_list))


# ---- 1. Read in rgSet
message(paste0('Reading input file from ', opt$rgset, '...\n'))

rgSet <- qread(opt$rgset)


methods <- unlist(strsplit(opt$methods, split=' '))
outputs <- unlist(strsplit(opt$outputs, split=' '))
reports <- unlist(strsplit(opt$reports, split=' '))

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