library(minfi, quietly=TRUE)
library(IlluminaHumanMethylationEPICmanifest, quietly=TRUE)
library(data.table, quietly=TRUE)
library(qs, quietly=TRUE)
library(optparse, quietly=TRUE)


# ---- 0. Parse CLI arguments
option_list <- list(
    make_option(c('-i', '--input'), 
        help=c("Path to the input RGChannelSet object."), 
        type="character"),
    make_option(c('-f', '--figures'),
        help='Path and filename to .pdf file qc figures will be written to.', 
        type='character'),
    make_option(c('-o', '--output'),
        help='Path and filename to save the per sample probe QC results to.', 
        type='character')
)

opt <- parse_args(OptionParser(option_list=option_list))


# ---- 1. Read in the RGSet from the previous step
message(paste0("Loading RGChannelSet object from ", opt$input, '...\n'))

rgSet <- qread(opt$input)


# ---- 2. Preprocess the RGChannelSet into a MethylSet by converting probe intensities to methylation values
message('Preprocessing RGChannelSet into MethylSet...\n')

methylSet <- preprocessRaw(rgSet)


# ---- 3. Run QC on 
message('Running initial QC on MethylSet and saving figures to ', opt$figures, '...\n')

methylSetQC <- getQC(methylSet)
pdf(opt$figures, width=6, height=6)
## TODO: Parameterize the bad sample cut-off
plotQC(methylSetQC, badSampleCutoff=9)
dev.off()


# ---- 4. Write bad samples to qc directory
meds <- (methylSetQC$mMed + methylSetQC$uMed)/2
badSampleIdxs <- which(meds < 9)
badSamples <- paste0("Bad sample idxs: ", paste0(badSampleIdxs, collapse=', '))
message(badSamples)

badSamplesDT <- data.table(
    'sample'=colnames(methylSet)[badSampleIdxs],
    'sample_idx'=badSampleIdxs
)
fwrite(badSamplesDT, file=paste0(opt$figures, '.badsamples.csv'))


## TODO: Add option to automatically drop bad samples

# ---- 4. Save methylSet
message('Saving MethylSet to ', opt$output, '...\n')

qsave(methylSet, file=opt$output)


message("Done!")