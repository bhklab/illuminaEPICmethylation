# ---- 0. Load dependencies
message('Loading script dependencies...\n')

# Suppress package load messages to ensure logs are not cluttered
suppressMessages({
    library(minfi, quietly=TRUE)
    library(data.table, quietly=TRUE)
    library(qs, quietly=TRUE)
    library(wateRmelon, quietly=TRUE)
})

# ---- 0. Parse Snakemake arguments
input <- snakemake@input
output <- snakemake@output

# ---- 1. Read in RGChannelSet
message(paste0('Reading in RGChannelSet from ', input$rgset, '...\n'))
rgSet <- qread(input$rgset)


# ---- 2. Get the bead count for each probe in each sample
message('Extracting bead count matrix...\n')
beadCountMatrix <- beadcount(rgSet)

message('Saving bead counts to: ', output$bead_counts, '...\n')
beadCountDT <- data.table(beadCountMatrix, keep.rownames='rownames')
fwrite(beadCountDT, file=output$bead_counts)

message("Filtering RGChannelSet to probes with more than three beads...\n")
perCount <- 1 / ncol(rgSet)

rgSetFiltered <- pfilter(mn=rgSet, un=rgSet, 
                    bc=beadCountMatrix, perCount=perCount,
                    perc=10.1, pthresh=100,
                    logical.return=TRUE)

print(rgSetFiltered)

message(paste0("Saving filtered RGChannelSet to: ", output$rgset_filtered, '...\n'))
qsave(rgSetFiltered, file=output$rgset_filtered)

message("Done!\n\n")
