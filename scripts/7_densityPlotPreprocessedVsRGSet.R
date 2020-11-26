# ---- 0. Load dependencies
message('Loading script dependencies...\n')

# Suppress package load messages to ensure logs are not cluttered
suppressMessages({
    library(minfi, quietly=TRUE)
    library(qs, quietly=TRUE)
    library(BiocParallel, quietly=TRUE)
    library(data.table, quietly=TRUE)
    library(ggplot2, quietly=TRUE)
    library(gridExtra, quietly=TRUE)
    library(grid, quietly=TRUE)
})

# ---- 0. Parse Snakemake arguments

input <- snakemake@input
output <- snakemake@output

methylSetPaths <- input$methylsets
plotTitles <- c('Unprocessed', input$plot_titles)
qcReportPath <- output$qc

# ---- 1. Read in methylSet
message("Reading in data from:\n\t", 
        paste0(input$rgset, '\n\t', paste0(methylSetPaths, collapse='\n\t'), '\n'))

methylSets <- bplapply(methylSetPaths, qread, nthreads=4)
rgSet <- qread(input$rgset, nthreads=20)

plotSets <- c(list(rgSet), methylSets)
rm(methylSets); rm(rgSet); gc()

# ---- 2. Define QC Functions

diffMethylVsUnmethylPeaks <- function(betaValues) {
    # Drop any NA values
    betaValues <- na.omit(betaValues)

    # Calculate the density values
    densityResults <- density(betaValues)
    densityBetas <- densityResults$x
    densityValues <- densityResults$y

    # Find the middle index
    middleIndex <- floor(length(densityValues) / 2)

    # Find the index of the max density above and below the middle index
    whichMaxMethylated <- middleIndex + 
        which.max(densityValues[seq(middleIndex + 1, length(densityValues))])
    whichMaxUnmethylated <- which.max(densityValues[seq_len(middleIndex)])

    # Get the peak of methylated vs unmethylated signal
    methylatedPeak <- densityBetas[whichMaxMethylated]
    unmethylatedPeak <- densityBetas[whichMaxUnmethylated]
    
    .countPeaks <- function(density) sum(diff(sign(diff(density))) == -2)
    
    numPeaks <- .countPeaks(densityValues)
    
    peakDifference <- round(methylatedPeak - unmethylatedPeak, 2)
    
    return(list('peakDiff'=peakDifference, 'numPeaks'=numPeaks))
}

calcMethylSignalOverNoise <- function(betaValues) {
    # Drop any NA values
    betaValues <- na.omit(betaValues)
    
    # Range of beta values is (0, 1), therefore AUC of density plot is 
    #   just number of obs over a specified interval
    methylAUC <- sum(betaValues >= 0.75)
    noiseAUC <- sum(betaValues < 0.75 & betaValues > 0.25)
    unmethylAUC <- sum(betaValues <= 0.25)

    signalOverNoise <- round((methylAUC + unmethylAUC) / noiseAUC, 2)
    return(signalOverNoise)
}

# ---- 3. Calculate QC Metrics and Build a DataTable

# Write functions to apply over each betaValue matrix
.sampleMethylVsUnmethylPeaks <- function(betaMatrix) {
    apply(betaMatrix, 2, FUN=diffMethylVsUnmethylPeaks)
}

.sampleMethylSignalOverNosie <- function(betaMatrix) {
    apply(betaMatrix, 2, FUN=calcMethylSignalOverNoise)
}

message("Extracting beta value matrices...\n")
# NOTE: This assumes the same samples for all plotSets
sampleNames <- make.names(colnames(plotSets[[1]]))


.getBetaNaOmit <- function(methylSet) na.omit(getBeta(methylSet))
betaMatrixList <- bplapply(plotSets, FUN=.getBetaNaOmit)
rm(plotSets); gc()

message("Calculating QC metrics...\n")
peakStats <- bplapply(betaMatrixList, FUN=.sampleMethylVsUnmethylPeaks)

.getPeakStat <- function(peakStats, statName) lapply(peakStats, `[[`, statName)

peakDifferences <- bplapply(peakStats, FUN=.getPeakStat, statName='peakDiff')
numPeaks <- bplapply(peakStats, FUN=.getPeakStat, statName='numPeaks')

aucRatios <- bplapply(betaMatrixList, FUN=.sampleMethylSignalOverNosie)

betaStatsDTs <- vector(mode='list', length=length(peakDifferences)) 
for (methodIdx in seq_along(peakDifferences)) {
    betaStatsDTs[[methodIdx]] <- data.table(
        'sample'=sampleNames,
        'peakDiff'=unlist(peakDifferences[[methodIdx]]),
        'numPeaks'=unlist(numPeaks[[methodIdx]]),
        'aucRatio'=unlist(aucRatios[[methodIdx]])
    )
}
names(betaStatsDTs) <- plotTitles
statsDT <- rbindlist(betaStatsDTs, idcol='preprocMethod')
fwrite(statsDT, file=qcReportPath)
rm(betaStatsDTs); gc()

# ---- 4. Assemble DataTable of Beta Values for Plotting

message('Building data.table object for plotting...\n')

betaDTList <- bplapply(betaMatrixList, as.data.table)
names(betaDTList) <- plotTitles
rm(betaMatrixList); gc()

.renameColumnsMakeNames <- function(DT) { 
    colnames(DT) <- make.names(colnames(DT)); return(DT)
}
betaDTList <- lapply(betaDTList, FUN=.renameColumnsMakeNames)

# ----- 5. Rendering Plots and Writing to PDF 
message(paste0("Writing plots to:\n\t", output$plots, '\n'))

.ggplotDensity <- function(DT, sampleName, plotTitle) {
    ggplotGrob(
        ggplot(data=DT, aes_string(sampleName)) +
            geom_density() + ggtitle(plotTitle) + xlab('Beta Values') +
            ylab("Density") + theme(plot.title=element_text(hjust=0.5))
    )
}

.densityPlotSample <- function(sampleName, betaDTList, plotTitles, statsDT) {
    message(paste0("Plotting: ", sampleName))
    
    grobs <- mapply(FUN=.ggplotDensity,
                    DT=betaDTList,
                    plotTitle=plotTitles,
                    MoreArgs=list('sampleName'=sampleName),
                    SIMPLIFY=FALSE)
    
    stats <- tableGrob(statsDT[sample == sampleName, -'sample'])
    
    allGrobs <- c(grobs, list(stats))
    return(allGrobs)
}

pdfGrobs <- lapply(sampleNames,
                     FUN=.densityPlotSample, 
                     betaDTList=betaDTList,
                     plotTitles=plotTitles,
                     statsDT=statsDT)

pdf(output$plots, height=11, width=8.5)

for (sampleIdx in seq_along(sampleNames)) {
    message(paste0('Writing ', sampleNames[sampleIdx], ' to pdf...'))
    title <- textGrob(sampleNames[sampleIdx], gp=gpar(fontsize=16, fontface='bold'))
    grid.arrange(grobs=pdfGrobs[[sampleIdx]], top=title)
}

dev.off()

message("Done!\n\n")
