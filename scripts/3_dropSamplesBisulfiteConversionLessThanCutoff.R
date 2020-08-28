# ---- 0. Load dependencies
message("Loading script dependencies...\n")

# Suppress package load messages to ensure logs are not cluttered
suppressMessages({
    library(minfi, quietly=TRUE)
    library(IlluminaHumanMethylationEPICmanifest, quietly=TRUE)
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
    make_option(c('-c', '--conversion_rate'),
        help='Bisulphite conversion cut-off for calling probes as detected, 
              as an integer specifying the percent converted.',
        type='numeric'),
    make_option(c('-P', '--path'),
        help='Path to save the output qc files to.',
        type='character'),
    make_option(c('-o', '--output'),
        help='Path and filename to save the filtered RGChannelSet to.', 
        type='character')
)

opt <- parse_args(OptionParser(option_list=option_list))


# ---- 1. Load RGChannelSet
message("Reading in RGChannelSet from: ", opt$rgset, '...\n')
rgSet <- qread(opt$rgset)


# ---- 2. Calculate the bisulphite conversion rates per sample and save to disk
message("Calculating sample bisulphite conversion rates...\n")
bisulphiteConversion <- bscon(rgSet)
bisulphiteConversionDT <- data.table(
                            'sample'=names(bisulphiteConversion),
                            'conversion_rate'=bisulphiteConversion
                            )


bsOut <- paste0(opt$path, '.bisulphite_conversions.csv')
fwrite(bisulphiteConversionDT, file=bsOut)


# ---- 3. Determine which how many and which samples fail bisulphite conversion qc
numPassed <- sum(bisulphiteConversionDT$conversion_rate > opt$conversion_rate)

message(paste0(numPassed, ' out of ', nrow(bisulphiteConversionDT), 
               ' samples passed qc at bisulphite conversion rate of ', 
               opt$conversion_rate, '%...\n'))

message(paste0('The following samples failed at bisulphite conversion rate of ', 
        opt$conversion_rate, '%:\n\t ', 
        paste0(bisulphiteConversionDT[conversion_rate < opt$conversion_rate]$sample, 
               collapse=',\n\t'), 
        '\n'))


# ---- 4. Drop samples failing the bisulphite conversion cut-off
message("Subsetting RGChannelSet to samples passing bisulphite QC...\n")
keepSamples <- bisulphiteConversion > opt$conversion_rate

rgSetFiltered <- rgSet[, keepSamples]
print(rgSetFiltered)
message('\n')

# ---- 5. Save RGChannelSet to disk
message(paste0("Saving RGChannel set to ", opt$output, '...\n'))
qsave(rgSetFiltered, file=opt$output)


message("Done!\n\n")
