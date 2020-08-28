library(minfi, quietly=TRUE)
library(optparse, quietly=TRUE)
library(data.table, quietly=TRUE)
library(BiocParallel, quietly=TRUE)
library(qs, quietly=TRUE)
library(plyr, quietly=TRUE)


# ---- 0. Parse CLI arguments
option_list <- list(
    make_option(c('-g', '--grsets'), 
        help='Path to write the file mapping from Beta value methylated region to CpG sites to.',
        type='character'),
    make_option(c('-d', '--methylation_data'),
        help='Path to read the per CpG methylation M and Beta values from',
        type='character'),
    make_option(c('-v', '--methylation_values'),
        help='Identifiers used in the paths for input files to separation of M and Beta values.',
        type='character'),
    make_option(c('-M', '--mappings_out'),
        help='Path to write the region to CpG mappings for M and Beta values to.',
        type='character')
    make_option(c('-V', '--methylation_out'),
        help='Path to write the methylated region M and Beta value to.',
        type='character')
)

opt <- parse_args(OptionParser(option_list=option_list))

.strsplitAndFlatten <- function(string, split) unlist(strsplit(string, split=split))

# Drop last value in opt since this is the 'help' flag
opt <- lapply(opt[-length(opt)],  # string argument 
              FUN=.strsplitAndFlatten, 
              split=' ')


# ---- 1. Read in GRSets and methylation data

grSets <- bplappy(opt$grsets, FUN=qread)
methyValues <- bplapply(opt$methylation_data, FUN=fread)


# ---- 2. Extract indexes per block

.getIndexes <- function(grSet) grSet[[2]]$indexes

indexes <- lapply(grSets, .getIndexes)


# ---- 3. Extract feature names

features <- lapply(methylValues, rownames)


# ---- 4. Map blocks to feature names (mappings)

.getFeaturesPerBlock <- function(block, features) features[block]

.mapFeaturesToRegions <- function(blocks, features) {
    lapply(blocks, FUN=.getFeaturesPerBlock, features=features)
}

mappings <- bplapply(indexes,  # blocks argument
                     FUN=.mapFeaturesToRegions, 
                     features=features)

.data.table.transpose <- function(x) data.table::transpose(data.table(x))

mappingDTs <- lapply(mappings, 

for (i in seq_along(mappings)) {
    fwrite(mappings[[i]], opt$mappings_out[i])
}


# ---- 5. Separate M and Beta value GRSets

isBeta <- grepl(pattern=sort(opt$methylation_values)[1], opt$grsets)

betaGRSets <- lapply(grSets[isBeta],  # x argument
                     FUN=`[[`,  # x[[i]] 
                     i=1)
mGRSets <- lapply(grSets[!isBeta], 
                  FUN=`[[`, 
                  i=1)

# ---- 6. Extract M and Beta Values

betaValues <- bplapply(betaGRSets, getBeta)
betaValueDTs <- bplapply(betaValues, 
                         FUN=data.table, 
                         keep.rownames='rownames')
mValues <- bplapply(mGRSets, getM)
mValueDTs <- bplapply(mValues, 
                      FUN=data.table, 
                      keep.rownames='rownames')