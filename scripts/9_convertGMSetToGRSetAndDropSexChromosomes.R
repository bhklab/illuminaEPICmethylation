# ---- 0. Load dependencies
message("Loading script dependencies...\n")

# Suppress package load messages to ensure logs are not cluttered
suppressMessages({
    library(minfi, quietly=TRUE)
    library(data.table, quietly=TRUE)
    library(qs, quietly=TRUE)
})

# ---- 0. Parse Snakemake arguments
input <- snakemake@input
output <- snakemake@output


# ---- 1. Read in data
message(paste0('Reading in MethylSet from: ', input$methylset, '...\n'))
methylSet <- qread(input$methylset)


# ---- 2. Annotate MethylSet to GenomicMethylSet and save 
message(paste0('Annotating MethylSet and saving GenomicMethylSet to: ', output$genomicmethylset, '...\n'))

gmSet <- mapToGenome(methylSet)
qsave(gmSet, output$genomicmethylset)


# ---- 3. Convert GMSet to GRSet
message(paste0('Converting GenomicMethylSet to GenomicRatioSet...\n'))
grSet <- ratioConvert(gmSet, what='both')


# ---- 4. Drop sex chromosomes and save to GenomicMemthylSet
message(paste0('Removing sex chromosomes from GenomicRatioSet and saving to: ', output$ratioset))

grSet_annotations <- getAnnotation(grSet)
keep <- !(grSet_annotations$chr %in% c('chrX', 'chrY'))
message('Keep features: ')
print(table(keep))

grSet <- grSet[keep, ]

qsave(grSet, output$ratioset)


message("Done!\n\n")
