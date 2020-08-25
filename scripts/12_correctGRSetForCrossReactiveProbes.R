library(minfi)
library(optparse)
library(maxprobes)
library(qs)


# ---- 0. Parse CLI arguments

option_list <- list(
    make_option(c('-g', '--grset'), 
        help=c("A comma separated list of paths to the directory for each microarray plate."),
        type="character"),
    make_option(c('-o', '--output'), 
        help='Path and filename to save the output to.', 
        type='character')
)

opt <- parse_args(OptionParser(option_list=option_list))


# ---- 1. Reading in GenomicRatioSet

grSet <- qread(opt$grset)


# ---- 2. Get list of EPIC cross-reactive probes
EPICxReactiveProbes <- maxprobes:::xreactive_probes()


# ---- 3. Filter out cross-reactive probes
keep <- !(featureNames(grSet) %in% EPICxReactiveProbes)
message("Keep: ")
print(table(keep))

grSetNoXReactive <- grSet[keep, ]


# ---- 4. Save filtered GenomicRatioSet
qsave(grSetNoXReactive, file=opt$output)