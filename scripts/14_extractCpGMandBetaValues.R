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

# ---- 1. Load in data
message(paste0("Reading in GenomicRatioSets from: \n\t", paste0(opt$grsets, collapse='\n\t'), '\n'))
grsets <- bplapply(opt$grsets, FUN=qread)

# ---- 2. Extract all beta values
message("Extracting Beta values...\n")
betaValues <- lapply(grsets, getBeta)
betaValueDTs <- lapply(betaValues, # object argument
                       FUN=data.table, 
                       keep.rownames='rownames')

# ---- 3. Save beta values as .csv
message(paste0("Saving Beta values to: \n\t", paste0(opt$beta_values, collapse='\n\t'), '\n'))
for (i in seq_along(betaValueDTs)) {
    fwrite(betaValueDTs[[i]], file=opt$beta_values[i])
}


# ---- 4. Extract all m values
message("Extracting M values...\n")
mValues <- lapply(grsets, getM)
mValueDTs <- lapply(mValues,  # object argument
                    FUN=data.table,
                    keep.rownames='rownames')


# ---- 5. Save m values as .csv
message(paste0("Saving M values to: \n\t", paste0(opt$beta_values, collapse=',\n\t'), '\n'))
for (i in seq_along(betaValueDTs)) {
    fwrite(mValueDTs[[i]], file=opt$m_values[i])
}


message("Done!\n\n")