# ---- 0. Load dependencies
message("Loading script dependencies...\n")

# Suppress package load messages to ensure logs are not cluttered
suppressMessages({
    library(minfi, quietly=TRUE)
    library(data.table, quietly=TRUE)
    library(BiocParallel, quietly=TRUE)
    library(qs, quietly=TRUE)
    library(plyr, quietly=TRUE)
})

# ---- 0. Parse CLI arguments
input <- snakemake@input
output <- snakemake@output

# ---- 1. Read in GRSets and methylation data
message(paste0("Reading in GRSets from:\n\t", paste0(input$grset, collapse='\n\t'), '\n'))
grSets <- bplapply(input$grsets, FUN=qread)


message(paste0("Reading in methylation values from:\n\t", 
               paste0(input$methylation_data, collapse='\n\t'), '\n'))
methylValues <- bplapply(input$methylation_data, FUN=fread)


# ---- 2. Extract indexes per block
message("Extracting methylated region indexes...\n")
.getIndexes <- function(grSet) grSet[[2]]$indexes

indexes <- lapply(grSets, .getIndexes)


# ---- 3. Extract feature names
message("Extracting feature names..\n")

.getDTrownames <- function(DT) DT$rownames
features <- lapply(methylValues, FUN=.getDTrownames)


# ---- 4. Map blocks to feature names (mappings)
message("Contructing feature to region mapping files...\n")
.getFeaturesPerBlock <- function(block, features) features[block]

.mapFeaturesToRegions <- function(blocks, features) {
    lapply(blocks, FUN=.getFeaturesPerBlock, features=features)
}

mappings <- bpmapply(FUN=.mapFeaturesToRegions,
                     blocks=indexes, features=features,  # FUN args
                     SIMPLIFY=FALSE)

mappingDTs <- lapply(mappings, FUN=data.table)                     

message(paste0("Writing methylated region mappings to:\n\t",
        paste0(output$mappings, collapse='\n\t'), '\n'))

for (i in seq_along(mappingDTs)) {
    fwrite(mappingDTs[[i]], output$mappings[i])
}


# ---- 5. Separate M and Beta value GRSets
message("Extracting M and Beta values from the GRSets...\n")
isBeta <- vapply(output$methyl_values, 
                 FUN=grepl, pattern='.*beta_values.*', 
                 FUN.VALUE=logical(1))

betaGRSets <- lapply(grSets[isBeta],  # x argument
                     FUN=`[[`,  # x[[i]] 
                     i=1)
mGRSets <- lapply(grSets[!isBeta], 
                  FUN=`[[`, 
                  i=1)

# ---- 6. Extract M and Beta Values

betaValues <- bplapply(betaGRSets, getBeta)
betaValueDTs <- bplapply(betaValues, 
                         FUN=data.table, 
                         keep.rownames=FALSE)

mValues <- bplapply(mGRSets, getM)
mValueDTs <- bplapply(mValues, 
                      FUN=data.table, 
                      keep.rownames=FALSE)


message(paste("Writing m and beta value data.tables to .csv: ", 
        paste0(output$methyl_values, collapse='\n\t'), '\n'))
for (i in seq_along(betaValueDTs)) {
    fwrite(betaValueDTs[[i]], file=output$methyl_values[isBeta][i])
    fwrite(mValueDTs[[i]], file=output$methyl_values[!isBeta][i])
}