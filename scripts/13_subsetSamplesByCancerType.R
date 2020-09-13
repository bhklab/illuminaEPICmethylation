# ---- 0. Load dependencies
message("Loading script dependencies...\n")

# Suppress package load messages to ensure logs are not cluttered
suppressMessages({
    library(minfi, quietly=TRUE)
    library(data.table, quietly=TRUE)
    library(qs, quietly=TRUE)
    library(optparse, quietly=TRUE)
    library(BiocParallel, quietly=TRUE)
})

# ---- 0. Parse CLI arguments
option_list <- list(
    make_option(c('-g', '--grset'), 
        help=c("A comma separated list of paths to the directory for each microarray plate."),
        type="character"),
    make_option(c('-s', '--subsets'), 
        help='Path and filename to save the output to.', 
        type='character'),
    make_option(c('-o', '--outputs'), 
        help='Path and filename to save the output to.', 
        type='character')
)

opt <- parse_args(OptionParser(option_list=option_list))


# ---- 1. Read in GenomicRatioSet
grSet <- qread(opt$grset)


# ---- 2. Match group names

# Parse string input into vector
cancerTypes <- unlist(strsplit(opt$subsets, split=' '))
cancerTypes <- lapply(cancerTypes, FUN=gsub, pattern='_', replacement=' ')

subsets <- lapply(cancerTypes,  # pattern
                  FUN=grepl,
                  x=colData(grSet)$Sample_Group)


# ---- 3. Subset GRSet based on 
.subsetCols <- function(j, x) x[i=TRUE, j]

grSets <- bplapply(subsets, 
                 FUN=.subsetCols, 
                 x=grSet)


# ---- 4. Save subsets to disk
message('Saving GenomicRatioSet subsets...\n')
outputs <- unlist(strsplit(opt$outputs, split=' '))
combinedGRSet <- outputs[1]
outputs <- outputs[-1]

qsave(grSet, file=combinedGRSet)

for (i in seq_along(grSets)) {
    message('...Saving :', outputs[i], '...')
    qsave(grSets[[i]], file=outputs[i])
}

message("Done!\n\n")