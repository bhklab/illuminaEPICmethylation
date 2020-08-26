library(minfi, quietly=TRUE)
library(optparse, quietly=TRUE)
library(BiocParallel, quietly=TRUE)
library(qs, quietly=TRUE)


# ---- 0. Parse CLI arguments
option_list <- list(
    make_option(c('-i', '--input'), 
        help='Path to normalized MethylSet.', 
        type='character'),
    make_option(c('-a', '--annotated'),
        help='Path to write the annotated GenomicMethylSet to.',
        type='character'),
    make_option(c('-o', '--output'),
        help='Path to write the GenomicRatioSet to.',
        type='character')
)

opt <- parse_args(OptionParser(option_list=option_list))


# ---- 1. Read in data
message(paste0('Reading in MethylSet from: ', opt$input, '...\n'))
methylSet <- qread(opt$input)


# ---- 2. Annotate MethylSet to GenomicMethylSet and save 
message(paste0('Annotating MethylSet and saving GenomicMethylSet to: ', opt$annotated, '...\n'))

gmSet <- mapToGenome(methylSet)
qsave(gmSet, opt$annotated)


# ---- 3. Convert GMSet to GRSet
message(paste0('Converting GenomicMethylSet to GenomicRatioSet...\n'))
grSet <- ratioConvert(gmSet, what='both')


# ---- 4. Drop sex chromosomes and save to GenomicMemthylSet
message(paste0('Removing sex chromosomes from GenomicRatioSet and saving to: ', opt$output))

grSet_annotations <- getAnnotation(grSet)
keep <- !(grSet_annotations$chr %in% c('chrX', 'chrY'))
message('Keep features: ')
print(table(keep))

grSet <- grSet[keep, ]

qsave(grSet, opt$output)


message("Done!\n\n")
