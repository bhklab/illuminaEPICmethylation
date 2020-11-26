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

# ---- 0. Parse Snakemake arguments
input <- snakemake@input
params <- snakemake@params
output <- snakemake@output

# ---- 1. Read in GenomicRatioSet
message("Reading in GenomicRatioSet from ", input$grset, "\n")
grSet <- qread(input$grset)


# ---- 2. Match group names
message("Performing subset by cancer type\n")

# Parse string input into vector
cancerTypes <- unlist(strsplit(params$subsets, split=' '))
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