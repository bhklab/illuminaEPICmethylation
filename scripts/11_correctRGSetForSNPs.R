library(minfi)
library(optparse)
library(BiocParallel)
library(qs)


# ---- 0. Parse CLI arguments

option_list <- list(
    make_option(c('-p', '--plates'), 
        help=c("A comma separated list of paths to the directory for each microarray plate."), 
        type="character"),
    make_option(c('-l', '--labels'), 
        help=c("A comma separated list of paths to the labels for each plate, in the same order as '--plates'."), 
        type="character"),
    make_option(c('-o', '--output'), 
        help='Path and filename to save the output to.', 
        type='character'),
    make_option(c('-n', '--nthread'), 
        help=c("The number of threads to parallelize over."), 
        type='integer', default=1)
)

opt <- parse_args(OptionParser(option_list=option_list))
