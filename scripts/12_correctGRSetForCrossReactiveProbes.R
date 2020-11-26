# ---- 0. Load dependencies
message("Loading script dependencies...\n")

# Suppress package load messages to ensure logs are not cluttered
suppressMessages({
    library(minfi, quietly=TRUE)
    library(data.table, quietly=TRUE)
    library(qs, quietly=TRUE)
    library(optparse, quietly=TRUE)
})

# ---- 0. Parse Snakemake arguments
input <- snakemake@input
output <- snakemake@output

# ---- 1. Reading in GenomicRatioSet
message("Reading in GenomcRatioSet from ", input$getset, '\n')
grSet <- qread(input$grset)


# ---- 2. Get list of EPIC cross-reactive probes
message("Filtering cross-reactive probes")
EPICxReactiveProbes <- maxprobes:::xreactive_probes()


# ---- 3. Filter out cross-reactive probes
keep <- !(featureNames(grSet) %in% EPICxReactiveProbes)
message("Keep: ")
print(table(keep))

grSetNoXReactive <- grSet[keep, ]

message("Saving cross-reactive probe filtered grSet to ", output$drop_xreactive_grset, '\n')
# ---- 4. Save filtered GenomicRatioSet
qsave(grSetNoXReactive, file=output$drop_xreactive_grset)

message("Done!\n\n")