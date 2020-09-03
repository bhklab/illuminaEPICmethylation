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
message('Parsing command line arguments...\n')
option_list <- list(
    make_option(c('-r', '--rgset'), 
        help="Path to the input RGChannelSet object.", 
        type="character"),
    make_option(c('-b', '--bead_counts'), 
        help="Path to save the bead count .csv file to.", 
        type="character"),
    make_option(c('-o', '--output'),
        help='Path and filename to save the filtered RGChannelSet to.', 
        type='character')
)

opt <- parse_args(OptionParser(option_list=option_list))


# ---- 1. Read in RGChannelSet
message(paste0('Reading in RGChannelSet from ', opt$input, '...\n'))
rgSet <- qread(opt$rgset)


# ---- 2. Get the bead count for each probe in each sample
message('Extracting bead count matrix...\n')
beadCountMatrix <- beadcount(rgSet)

message('Saving bead counts to: ', opt$bead_counts, '...\n')
beadCountDT <- data.table(beadCountMatrix, keep.rownames='rownames')
fwrite(beadCountDT, file=opt$bead_counts)

message("Filtering RGChannelSet to probes with more than three beads...\n")
perCount <- 1 / ncol(rgSet)

rgSetFiltered <- pfilter(mn=rgSet, un=rgSet, 
                    bc=beadCountMatrix, perCount=perCount,
                    perc=10.1, pthresh=100,
                    logical.return=TRUE)

print(rgSetFiltered)

message(paste0("Saving filtered RGChannelSet to: ", opt$output, '...\n'))
qsave(rgSetFiltered, file=opt$output)


message("Done!\n\n")
