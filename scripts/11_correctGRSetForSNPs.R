# ---- 0. Load dependencies
message("Loading script dependencies...\n")

# Suppress package load messages to ensure logs are not cluttered
suppressMessages({
    library(minfi, quietly=TRUE)
    library(data.table, quietly=TRUE)
    library(qs, quietly=TRUE)
    library(BiocParallel, quietly=TRUE)
})

# ---- 0. Parse CLI arguments
input <- snakemake@input
output <- snakemake@output

# ---- 1. Load filtered GenomicRatioSet
message(paste0("Loading filtered GenomicRatioSet from: ", input$grset, '...\n'))

grSet <- qread(input$grset)


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
message(paste0("Writing SNP filtered GenomicRatioSet to: ", output$drop_snps_grset, '...\n'))

qsave(grSetDropSNPs, file=output$drop_snps_grset)


message("Done!\n\n")