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


# ---- 0. Parse Snakemake arguments

input <- snakemake@input
params <- snakemake@params
output <- snakemake@output

# ---- 1. Load RGChannelSet
message("Reading in RGChannelSet from: ", input$rgset, '...\n')
rgSet <- qread(input$rgset)


# ---- 2. Calculate the bisulphite conversion rates per sample and save to disk
message("Calculating sample bisulphite conversion rates...\n")
bisulphiteConversion <- bscon(rgSet)
bisulphiteConversionDT <- data.table(
                            'sample'=names(bisulphiteConversion),
                            'conversion_rate'=bisulphiteConversion
                            )


fwrite(bisulphiteConversionDT, file=output$bisulphite_qc)


# ---- 3. Determine which how many and which samples fail bisulphite conversion qc
conversion_minimum <- params$bisulphite_conversion_rate

numPassed <- sum(bisulphiteConversionDT$conversion_rate > conversion_minimum)

message(paste0(numPassed, ' out of ', nrow(bisulphiteConversionDT), 
               ' samples passed qc at bisulphite conversion rate of ', 
               conversion_minimum, '%...\n'))

message(paste0('The following samples failed at bisulphite conversion rate of ', 
        conversion_minimum, '%:\n\t ', 
        paste0(bisulphiteConversionDT[conversion_rate < conversion_minimum]$sample, 
               collapse=',\n\t'), 
        '\n'))


# ---- 4. Drop samples failing the bisulphite conversion cut-off
message("Subsetting RGChannelSet to samples passing bisulphite QC...\n")
keepSamples <- bisulphiteConversion > coversion_minimum

rgSetFiltered <- rgSet[, keepSamples]
print(rgSetFiltered)
message('\n')

# ---- 5. Save RGChannelSet to disk
message(paste0("Saving RGChannel set to ", output$rgset_filtered, '...\n'))
qsave(rgSetFiltered, file=output$rgset_filtered)


message("Done!\n\n")
