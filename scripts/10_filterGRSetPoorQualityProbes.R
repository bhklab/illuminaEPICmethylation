library(minfi)
library(optparse)
library(data.table)
library(matrixStats)
library(qs)


# ---- 0. Parse CLI arguments
option_list <- list(
    make_option(c('-g', '--grset'), 
        help='Path to GenomicRatioSet.', 
        type='character'),
    make_option(c('-p', '--pvalues'), 
        help='Path to QC3 detection p-value .csv file.', 
        type='character'),
    make_option(c('-o', '--output'),
        help='Path to write the filtered GenomicRatioSet to.',
        type='character')
)

opt <- parse_args(OptionParser(option_list=option_list))

# ---- 1. Load data
message(paste0('Loading RGSet from: ', opt$input, '...\n'))

grSet <- qread(opt$grset)
pValues <- fread(opt$pvalues)


# --- 2 . Filter poor quality probes
message(paste0('Filtering probes failing...\n'))

# Match probes in pValues to rgSet
pValues <- pValues[rownames %in% featureNames(grSet), ]
keep <- rowSums(pValues[, -'rownames'] < 0.01) == ncol(grSet)
message("Keep: \n")
print(table(keep))

grSet <- grSet[keep, ]


# ---- 3. Save filtered RGSet
message(paste0('\nSaving GenomicRatioSet to: ', opt$output, '...\n'))

qsave(grSet, opt$output)

message("Done...\n\n")