# ---- 0. Load dependencies
message("Loading script dependencies...\n")

# Suppress package load messages to ensure logs are not cluttered
suppressMessages({
    library(minfi, quietly=TRUE)
    library(data.table, quietly=TRUE)
    library(qs, quietly=TRUE)
    library(optparse, quietly=TRUE)
    library(wateRmelon, quietly=TRUE)
})


# ---- 0. Parse CLI arguments
message("Parsing command line arguments...\n")
option_list <- list(
    make_option(c('-r', '--rgset'), 
        help="Path to the input RGChannelSet object.", 
        type="character"),
    make_option(c('-o', '--output'),
        help='Path and filename to save the filtered RGChannelSet to.', 
        type='character')
)

opt <- parse_args(OptionParser(option_list=option_list))


# ---- 1. Read in RGChannelSet

rgSet <- qread(opt$rgset)


# ---- 2. Get the bead count for each probe in each sample

beadCountMatrix <- beadcount(rgSet)
