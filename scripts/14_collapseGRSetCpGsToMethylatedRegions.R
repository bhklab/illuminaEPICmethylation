library(minfi, quietly=TRUE)
library(optparse, quietly=TRUE)
library(data.table, quietly=TRUE)
library(BiocParallel, quietly=TRUE)
library(qs, quietly=TRUE)

# ---- 0. Parse CLI arguments
option_list <- list(
    make_option(c('-g', '--grsets'), 
        help=c("A comma separated list of paths to the directory for each microarray plate."),
        type="character"),
    make_option(c('-b', '--beta_region_grsets'), 
        help='Path to write the file mapping from Beta value methylated region to CpG sites to.',
        type='character'),
    make_option(c('-m', '--m_region_grsets'),
        help='Path to write the file mapping from Beta value methylated region to CpG sites to.',
        type='character')
)

opt <- parse_args(OptionParser(option_list=option_list))

.strsplitAndFlatten <- function(string, split) unlist(strsplit(string, split=split))

# Drop last value in opt since this is the 'help' flag
opt <- lapply(opt[-length(opt)],  # string argument 
              FUN=.strsplitAndFlatten, 
              split=' ')


# ---- 1. Read in grsets
message(paste0("Reading in GenomicRatioSets from: \n\t", paste0(opt$grsets, collapse=',\n\t'), '\n'))
grSets <- bplapply(opt$grsets, FUN=qread)


# ---- 2. Collapse adjacent CpG sites for Beta values
message("Collapsing adjacent CpG sites for Beta values...\n")
grSetsBetaRegions <- bplapply(grSets,  # object argument
                              FUN=cpgCollapse,
                              what='Beta')


# ---- 3. Collapse adjacent Cpg sites for M values
message("Collapsing adjacent CpG sites for M values...\n")
grSetsMRegions <- bplapply(grSets,  # object arugment
                           FUN=cpgCollapse,
                           what='M')


# ---- 4. Save collapsed GenomicRatioSets
message(paste0("Saving collapsed GenomicRatioSets for Beta values to: ", paste0(opt$beta_region_grsets, collapse=',\n\t'), '\n'))
message(paste0("Saving collapsed GenomicRatioSets for M values to: ", paste0(opt$m_region_grsets, collapse=',\n\t'), '\n'))
for (i in seq_along(grSetsBetaRegions)) {
    qsave(grSetsBetaRegions[[i]], file=opt$beta_region_grsets[i])
    qsave(grSetsMRegions[[i]], file=opt$m_region_grsets[i])
}