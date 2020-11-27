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

# ---- 1. Read in grsets
message(paste0("Reading in GenomicRatioSets from: \n\t", paste0(input$grsets, collapse=',\n\t'), '\n'))
grSets <- bplapply(input$grsets, FUN=qread)


# ---- 2. Collapse adjacent CpG sites for Beta values
message("Collapsing adjacent CpG sites for Beta values...\n")
grSetsBetaRegions <- bplapply(grSets,  # object argument
                              FUN=cpgCollapse,
                              what='Beta')


# ---- 3. Collapse adjacent Cpg sites for M values
message("Collapsing adjacent CpG sites for M values...\n")
grSetsMRegions <- bplapply(grSets,  # object arugment
                           FUN=cpgCollapse,
                           what='M')


# ---- 4. Save collapsed GenomicRatioSets
message(paste0("Saving collapsed GenomicRatioSets for Beta values to: ", paste0(output$beta_region_grsets, collapse=',\n\t'), '\n'))
message(paste0("Saving collapsed GenomicRatioSets for M values to: ", paste0(output$m_region_grsets, collapse=',\n\t'), '\n'))
for (i in seq_along(grSetsBetaRegions)) {
    qsave(grSetsBetaRegions[[i]], file=output$beta_region_grsets[i])
    qsave(grSetsMRegions[[i]], file=output$m_region_grsets[i])
}