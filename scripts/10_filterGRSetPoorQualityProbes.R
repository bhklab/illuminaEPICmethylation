# ---- 0. Load dependencies
message("Loading script dependencies...\n")

# Suppress package load messages to ensure logs are not cluttered
suppressMessages({
    library(minfi, quietly=TRUE)
    library(data.table, quietly=TRUE)
    library(qs, quietly=TRUE)
    library(matrixStats, quietly=TRUE)
})

# ---- 0. Parse CLI arguments
input <- snakemake@input
output <- snakemake@output


# ---- 1. Load data
message(paste0('Loading RGSet from: ', input$grset, '...\n'))

grSet <- qread(input$grset)
pValues <- fread(input$pvalues)

# --- 2 . Filter poor quality probes
message(paste0('Filtering probes failing...\n'))

# Match probes in pValues to rgSet
colnames(grSet) <- trimws(colnames(grSet))
pValuesSub <- pValues[rownames %in% rownames(grSet), .SD, .SDcol=colnames(pValues) %in% c('rownames', colnames(grSet))]
keep <- rowSums(pValuesSub[, -'rownames'] < 0.01) == ncol(grSet)
message("Keep: \n")
print(table(keep))

grSet <- grSet[keep, ]


# ---- 3. Save filtered RGSet
message(paste0('\nSaving GenomicRatioSet to: ', output$filtered_grset, '...\n'))

qsave(grSet, output$filtered_grset)


message("Done...\n\n")