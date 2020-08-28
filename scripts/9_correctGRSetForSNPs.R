library(minfi, quietly=TRUE)
library(optparse, quietly=TRUE)
library(BiocParallel, quietly=TRUE)
library(qs, quietly=TRUE)


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


# ---- 1. Load filtered GenomicRatioSet
message(paste0("Loading filtered GenomicRatioSet from: ", opt$grset, '...\n'))

grSet <- qread(opt$grset)


# ---- 2. Extract SNPs
message("Extracting SNP information...\n")

SNPs <- getSnpInfo(grSet)

message("Number of CpG SNPs: ")
print(sum(grepl(SNPs$CpG_rs, pattern='rs*')))

message("Number of SBE SNPs: ")
print(sum(grepl(SNPs$SBE_rs, pattern='rs*')))


# ---- 3. Adding SNP info to GenomicRatioSet
message("Adding SNP info to GenomicRatioSet object...\n")

grSetSNPs <- addSnpInfo(grSet)


# ---- 4. Drop probes containing SNP at CpG interrogation or at single NT extension
message("Dropping loci containing SNPs...\n")
grSetDropSNPs <- dropLociWithSnps(grSetSNPs, snps=c('SBE', 'CpG'), maf=0)


# ---- 5. Writing SNP filtered GenomicRatioSet
message(paste0("Writing SNP filtered GenomicRatioSet to: ", opt$output, '...\n'))

qsave(grSetDropSNPs, file=opt$output)


message("Done!\n\n")