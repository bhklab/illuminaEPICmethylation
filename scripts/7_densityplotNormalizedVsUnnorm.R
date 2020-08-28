library(minfi, quietly=TRUE)
library(optparse, quietly=TRUE)
library(RColorBrewer, quietly=TRUE)
library(qs, quietly=TRUE)


# ---- 0. Parse CLI arguments
option_list <- list(
    make_option(c('-q', '--qualityControl'), 
        help=c("Paths to the files containing the RGChannelSets from the two quality 
               control steps preceeding normalization."), 
        type="character"),
    make_option(c('-n', '--normalized'), 
        help=c("Path to the GenomicMethylSet containging the normalzied output from QC step 3."), 
        type='character'),
    make_option(c('-o', '--output'), 
        help='Path to write the pdf file containing the normalized vs unnormalized QC2 plots',
        type='character')
)

opt <- parse_args(OptionParser(option_list=option_list))

qc <- unlist(strsplit(opt$qualityControl, split=' '))
output <- unlist(strsplit(opt$output, split= ' '))

if (length(qc) != length(output)) stop("Must include one input qc file per output plot!")


# ---- 1. Read in methylation data from QC2, QC3 and normalzied QC3
message(paste0("Loading data for QC plots from: ", 
        paste0(c(opt$plates, opt$qualityControl2, opt$qualityControl3), 
                 collapse=', '), '...\n'))

rgSetQc2 <- qread(qc[[1]])
rgSetQc3 <- qread(qc[[2]])
methylSet <- qread(opt$normalized)


# ---- 3. Write the QC2 vs Normalized Plots to Pdf 
message(paste0("Writing normalized vs unnormalzied QC2 plots to : ", opt$output1, '...\n'))

pdf(output[[1]])
par(mfrow=c(1,2))

densityPlot(rgSetQc2, sampGroups=colData(rgSetQc2)$Sample_Group, 
    main='QC2', legend=FALSE, pal=brewer.pal(9, "Set1"))
legend('topleft', legend=levels(factor(colData(rgSetQc2)$Sample_Group)),
      text.col=brewer.pal(9, "Set1"))

densityPlot(getBeta(methylSet), sampGroups=colData(rgSetQc2)$Sample_Group, 
    main='QC2 Normalized', legend=FALSE, pal=brewer.pal(9, "Set1"))
legend('topleft', legend=levels(factor(colData(rgSetQc2)$Sample_Group)),
      text.col=brewer.pal(9, "Set1"))

dev.off()


# ---- 4. Write the QC3 vs Noramalized Plots to Pdf
message(paste0("Writing normalized vs unnormalzied QC3 plots to : ", opt$output2, '...\n'))

pdf(output[[2]])
par(mfrow=c(1,2))

densityPlot(rgSetQc3, sampGroups=colData(rgSetQc3)$Sample_Group, 
    main='QC3', legend=FALSE, pal=brewer.pal(8, "Set1"))
legend('topleft', 
       legend=levels(factor(colData(rgSetQc3)$Sample_Group)),
       text.col=brewer.pal(8, "Set1"))

densityPlot(getBeta(methylSet), sampGroups=colData(rgSetQc3)$Sample_Group, 
    main='QC3 Normalized', legend=FALSE, pal=brewer.pal(8, "Set1"))
legend('topleft', 
       legend=levels(factor(colData(rgSetQc3)$Sample_Group)),
       text.col=brewer.pal(8, "Set1"))

dev.off()


message("Done!\n\n")