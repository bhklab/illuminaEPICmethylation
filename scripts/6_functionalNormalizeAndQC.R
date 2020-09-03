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
    make_option(c('-i', '--input'),
        help=c("Path to the quality controlled RGChannelSet object."),
        type="character"),
    make_option(c('-o', '--output'),
        help=c("Path to save the resulting methylSet to."),
        type="character"),
    make_option(c('-r', '--report'), 
        help='Path to save the qc file to.', 
        type='character')
)

opt <- parse_args(OptionParser(option_list=option_list))


# ---- 1. Read in rgSet
message(paste0('Reading input file from ', opt$input, '...\n'))

rgSet <- qread(opt$input)


# ---- 2. Functional normalize the rgSet to a methylSet
message(paste0('Functional normalzing the rgSet...\n'))

methylSet <- preprocessFunnorm(rgSet, ratioConvert=FALSE)
rm(rgSet); gc()

# ---- 3. Run QC on the methylSet
message('Running minfiQC on the methylSet...\n')

qcResults <- minfiQC(methylSet, verbose=TRUE)


# ---- 4. Extract the resulting methylSet and QC file
message("Extracting QC results...\n")

methylSet <- qcResults[[1]]

qcDf <- qcResults[[2]]


# ---- 5. Write output files
message(paste0("Saving methylSet to ", opt$output, '...\n'))

qsave(methylSet, file=opt$output)

message(paste0("Saving QC results to .csv in ", opt$report, '...\n'))
write.csv(qcDf, file=opt$report)


message("Done!\n\n")