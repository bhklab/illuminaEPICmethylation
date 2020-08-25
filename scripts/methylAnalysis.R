###
# Illumina Epic LMS Methylation Analysis with minfi
# by Christopher Eeles (CE) & Joanna Pryzbyl (JP)
#
# Adapted from:
# "A cross-package Bioconductor workflow for analysing methylation array data" from Bioconductor
# - Available at: https://www.bioconductor.org/packages/devel/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html
# "BioC2014 - Minfi tutorial"
# - Available at: https://www.bioconductor.org/help/course-materials/2014/BioC2014/minfi_BioC2014.pdf
#
###

###########
# Current unknowns/challenges:
#
# - Determine how to adjust for DNA degradation from FFPE samples
#   - Has this already been accounted for in the experimental pipeline?
#   - What further steps should be taken?
#
# - What is the best per feature metric to compare in Similarity Network Fusion analysis?
#   - Current assumption is to use genomic regions, as they may be more biologically interpretable
#
# - Sex chromosomes are usually the largest source of variation in methylation analysis
#   - Current assumption is to remove sex chromosomes
#
###########

library(minfi)
library(minfiData)
library(maxprobes)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICmanifest)
library(data.table)
library(RColorBrewer)
library(matrixStats)
library(qs)

# Clear any objects in environment to save on RAM
rm(list=ls())

####
# Data Import
####

print("Configuring directories...")

# Configure path to data
data_dir <- file.path('..', 'data', 'raw_data')
metadata_dir <- file.path('..', 'data', 'metadata')
results_dir <- file.path('..', 'analysis', '1_pipeline_results')
qc_dir <- file.path('..', 'analysis', '2_pipeline_qc')

# Get nthreads for parallelizing saves
nthread <- parallel::detectCores()

## Read data by plate

print("Reading data...")

# Configure plate paths
p1_dir <- file.path(data_dir, "Plate1")
p2_dir <- file.path(data_dir, "Plate2")
p3_dir <- file.path(data_dir, "Plate3")

# Create methylation array sheets
p1_sheet <- read.metharray.sheet(p1_dir)
p2_sheet <- read.metharray.sheet(p2_dir)
p3_sheet <- read.metharray.sheet(p3_dir)

# Import plate labels
## NOTE: Strings as factors should be FALSE by default in R 4.0
p1_labels <- read.csv(file.path(metadata_dir, "labels", "plate1_labels.csv"), stringsAsFactors = FALSE)
p2_labels <- read.csv(file.path(metadata_dir, "labels", "plate2_labels.csv"), stringsAsFactors = FALSE)
p3_labels <- read.csv(file.path(metadata_dir, "labels", "plate3_labels.csv"), stringsAsFactors = FALSE)

# Label sample groups
p1_sheet$Sample_Group <- p1_labels$Group
p2_sheet$Sample_Group <- p2_labels$Group
p3_sheet$Sample_Group <- p3_labels$Group

# Read data into rgSets
p1_RG <- read.metharray.exp(targets = p1_sheet)
p2_RG <- read.metharray.exp(targets = p2_sheet)
p3_RG <- read.metharray.exp(targets = p3_sheet)

# Label RGSet objects using methylation array sheets
sampleNames(p1_RG) <- paste(p1_sheet$Sample_Plate, p1_sheet$Sample_Well,
  p1_sheet$Sample_Group, sep = "-")
sampleNames(p2_RG) <- paste(p2_sheet$Sample_Plate, p2_sheet$Sample_Well,
  p2_sheet$Sample_Group, sep = "-")
sampleNames(p3_RG) <- paste(p3_sheet$Sample_Plate, p3_sheet$Sample_Well,
  p3_sheet$Sample_Group, sep = "-")

## Add metadata
p1_RG$metadata <- p1_sheet
p2_RG$metadata <- p2_sheet
p3_RG$metadata <- p3_sheet

print("Merging data...")


## Merge the plates into one rgSet
p12_RG <- combineArrays(p1_RG, p2_RG, outType = "IlluminaHumanMethylationEPIC")
stanford_meth_rgSet <- combineArrays(p12_RG, p3_RG, outType = "IlluminaHumanMethylationEPIC")

## Save stanford_meth_rgSet object
qsave(stanford_meth_rgSet, file = file.path(results_dir, "1_stanford_meth_rgSet_raw.qs"), nthread=nthread)
rm(p1_RG, p2_RG, p3_RG, p12_RG)
gc()

####
# QC
####

print("Starting QC...")

## Calculate p values for all probes for all samples
stanford_meth_pvals <- detectionP(stanford_meth_rgSet)

class(stanford_meth_pvals)
# [1] "matrix"

dim(stanford_meth_pvals)
# [1] 866091    304

# Identify genomic positions with p value > 0.01 ("failed" genomic positions)
failed <- stanford_meth_pvals > 0.01

# Calculate fraction of failed positions per sample
fractions_of_failed_positions <- colMeans(failed)

# Calculate a sum of failed positions per sample
sums_of_failed_positions <- colSums(failed)

# alternative way to calculate this:
# apply(stanford_meth_pvals, 2, sum)/nrow(stanford_meth_pvals)

# Save a file with the fraction of failed positions per sample
write.table(fractions_of_failed_positions, 
            file=file.path(qc_dir, "pval01_fractions_of_failed_positions_per_sample_304samples.tsv"), 
            sep = "\t ", na = "NA", row.names = TRUE, col.names = TRUE)

# Save a file with the fraction of failed positions per sample
write.table(sums_of_failed_positions , 
            file=file.path(qc_dir, "pval01_sums_of_failed_positions_per_sample_304samples.tsv"),
            sep = "\t ", na = "NA", row.names = TRUE, col.names = TRUE)

# How many genomic positions failed  p value QC in >50% / >10% / >5% / >1% of samples? At least 1 sample?
sum(rowMeans(failed) > 0.5) 
# [1] 589
sum(rowMeans(failed) > 0.1) 
# [1] 18579
sum(rowMeans(failed) > 0.05) 
# [1] 62729
sum(rowMeans(failed) > 0.01) 
# [1] 546068
sum(rowMeans(failed) > 0.0033) 
# [1] 801597

print("...QC report")

## Generate QC Report on all 304 samples
qcReport(
  stanford_meth_rgSet,
  sampNames = sampleNames(stanford_meth_rgSet),
  sampGroups = colData(stanford_meth_rgSet)$Sample_Group,
  pdf = file.path(qc_dir, "stanford_meth_QC1_report.pdf")
)

## JP: Interpretation of this qcReport file:
## The densityPlot function produces density plots of the methylation Beta values for all
## samples, typically colored by sample group. While the density plots in Figure 1 are useful
## for identifying deviant samples, it is not easy to identify the specific problem sample. If
## there is a concern about outlier samples, a useful follow-up is the “bean” plot (Figure 2)
## that shows each sample in its own section. While the shape of the distribution for “good”
## samples will differ from experiment to experiment, many conditions have methylation profiles
## characterized by two modes - one with close to 0% methylation, and a second at close to
## 100% methylation.


## Another QC plot recommended by the authors of minfi package
## minfi provides a simple quality control plot that uses the log median intensity in both the
## methylated (M) and unmethylated (U) channels. When plotting these two medians against
## each other, it has been observed that good samples cluster together, while failed samples
## tend to separate and have lower median intensities [1]. In general, we advice users to make
## the plot and make a judgement. The line separating ”bad” from ”good” samples represent a
## useful cutoff, which may have to be adapted to a specific dataset. 

## We need to apply preprocessRaw() first to the dataset:
## In order to obtain the methylated and unmethylated signals, we need to convert the RGChannelSet to an object containing 
## the methylated and unmethylated signals using the function preprocessRaw. 
## It takes as input a RGChannelSet and converts the red and green intensities to methylated and unmethylated signals 
## according to the special 450K probe design, and returns the converted signals in a new object of class MethylSet. 
## It does not perform any normalization.

MSet <- preprocessRaw(stanford_meth_rgSet)
MSet
qc <- getQC(MSet)
pdf(file.path(qc_dir, "preprocessed_stanford_meth_rgSet_qcPlot.pdf"))
plotQC(qc, badSampleCutoff = 9) # some sample are outsie the range of this plot
dev.off()

## shinyMethyl is particularly useful to visualize all plots at the same time in an interactive fashion.
## FIXME:: Can't install shiny on H4H cluster?
# myShinyMethylSet = shinySummarize(stanford_meth_rgSet)
# qsave(file.path(qc_dir, "myShinyMethylSet.qs"), nthread=nthread)

print("...QC1")

## Outlier samples with low signal from methylated and unmethylated channels, low values of control probes and poor density plots:
## Remove 15 samples that failed at least two different QC metrics

p1_failed_expQC <- c("B5", "C2", "D2", "F6", "H1", "H3")
p2_failed_expQC <- c("D5", "D6", "F5", "F6", "G6", "H5")
p3_failed_expQC <- c("C7", "C8", "D5")
p4_failed_expQC <- c()

stanford_meth_rgSet_QC1 <- stanford_meth_rgSet[,
 which(!(
   (colData(stanford_meth_rgSet)$Sample_Plate == "Plt-1" &
     (colData(stanford_meth_rgSet)$Sample_Well %in% p1_failed_expQC)) |
   (colData(stanford_meth_rgSet)$Sample_Plate == "Plt2" &
     (colData(stanford_meth_rgSet)$Sample_Well %in% p2_failed_expQC)) |
(colData(stanford_meth_rgSet)$Sample_Plate == "Plt3" &
  (colData(stanford_meth_rgSet)$Sample_Well %in% p3_failed_expQC)) # |
# (colData(stanford_meth_rgSet)$Sample_Plate == "Plt4" &
#   (colData(stanford_meth_rgSet)$Sample_Well %in% p4_failed_expQC))
 ))
]

stanford_meth_rgSet_QC1
# class: RGChannelSet 
# dim: 1051815 289 
# metadata(0):
#   assays(2): Green Red
# rownames(1051815): 1600101 1600111 ... 99810990 99810992
# rowData names(0):
#   colnames(289): Plt-1-A1-Control Plt-1-B1-DTF_het ... Plt4-G2-GIST Plt4-H2-LM_old
# colData names(11): Sample_Name Sample_Well ... metadata ArrayTypes
# Annotation
# array: IlluminaHumanMethylationEPIC
# annotation: ilm10b4.hg19

# Save the dataset with the outlier samples removed
qsave(stanford_meth_rgSet_QC1, file = file.path(results_dir, "2_stanford_meth_rgSet_QC1.qs"), nthread=nthread)

print("...QC2")

## Check the numbers of failed p values in samples that are replicate and remove the replicate that has higher number of failed p values
p1_failed_expQC <- c("B1", "C8")
p2_failed_expQC <- c("A6", "E7", "E11")
p3_failed_expQC <- c()
p4_failed_expQC <- c()

stanford_meth_rgSet_QC2 <- stanford_meth_rgSet_QC1[,
                                               which(!(
                                                 (colData(stanford_meth_rgSet_QC1)$Sample_Plate == "Plt-1" &
                                                    (colData(stanford_meth_rgSet_QC1)$Sample_Well %in% p1_failed_expQC)) |
                                                   (colData(stanford_meth_rgSet_QC1)$Sample_Plate == "Plt2" &
                                                      (colData(stanford_meth_rgSet_QC1)$Sample_Well %in% p2_failed_expQC)) # |
                                                #  (colData(stanford_meth_rgSet)$Sample_Plate == "Plt3" &
                                                #  (colData(stanford_meth_rgSet)$Sample_Well %in% p3_failed_expQC)) |
                                                 # (colData(stanford_meth_rgSet)$Sample_Plate == "Plt4" &
                                                 #   (colData(stanford_meth_rgSet)$Sample_Well %in% p4_failed_expQC))
                                               ))
                                               ]

stanford_meth_rgSet_QC2
# class: RGChannelSet
# dim: 1051815 284 
# metadata(0):
#   assays(2): Green Red
# rownames(1051815): 1600101 1600111 ... 99810990 99810992
# rowData names(0):
#   colnames(284): Plt-1-A1-Control Plt-1-C1-GIST ... Plt4-G2-GIST Plt4-H2-LM_old
# colData names(11): Sample_Name Sample_Well ... metadata ArrayTypes
# Annotation
# array: IlluminaHumanMethylationEPIC
# annotation: ilm10b4.hg19

# Save the dataset with the outlier samples removed
qsave(stanford_meth_rgSet_QC2, file = file.path(results_dir, "2_stanford_meth_rgSet_QC2.qs"), nthread=nthread)

print("...QC3")

## Calculate p values for all probes for the remaining samples
stanford_meth_pvals <- detectionP(stanford_meth_rgSet_QC2)

class(stanford_meth_pvals)
# [1] "matrix"

dim(stanford_meth_pvals)
# [1] 866091    284

# Identify genomic positions with p value > 0.01 ("failed" genomic positions)
failed <- stanford_meth_pvals > 0.01

# Calculate fraction of failed positions per sample
fractions_of_failed_positions <- colMeans(failed)

# Calculate a sum of failed positions per sample
sums_of_failed_positions <- colSums(failed)

# alternative way to calculate this:
# apply(stanford_meth_pvals, 2, sum)/nrow(stanford_meth_pvals)

# Save a file with the fraction of failed positions per sample
write.table(fractions_of_failed_positions,
            file=file.path(qc_dir, "pval01_fractions_of_failed_positions_per_sample_284samples.tsv"),
            sep="\t ", na="NA", row.names=TRUE, col.names=TRUE)

# Save a file with the fraction of failed positions per sample
write.table(sums_of_failed_positions, 
            file=file.path(qc_dir, "pval01_sums_of_failed_positions_per_sample_284samples.tsv"),
            sep = "\t ", na = "NA", row.names = TRUE, col.names = TRUE)

# How many genomic positions failed  p value QC in >50% / >10% / >5% / >1% of samples? At least 1 sample (1/284 = 0.00352113)?
sum(rowMeans(failed) > 0.5) 
# [1] 559
sum(rowMeans(failed) > 0.1) 
# [1] 11170
sum(rowMeans(failed) > 0.05) 
# [1] 28404
sum(rowMeans(failed) > 0.01) 
# [1] 137297
sum(rowMeans(failed) > 0.0036) 
# [1] 180675


## Remove 22 Jurkat DNA samples (these are just controls for assay performance)
## I think fully methylated Jurkat DNA may affect normalization? Agreed
p1_failed_expQC <- c("B2", "D5", "F9", "A1", "C4", "E8", "G11")
p2_failed_expQC <- c("F3", "E6", "C8", "A12", "E4", "D7", "B11")
p3_failed_expQC <- c("B7", "H2", "G3", "D5", "A6", "B9", "A10")
p4_failed_expQC <- c("A1")

stanford_meth_rgSet_QC3 <- stanford_meth_rgSet_QC2[,
                                                   which(!(
                                                     (colData(stanford_meth_rgSet_QC2)$Sample_Plate == "Plt-1" &
                                                        (colData(stanford_meth_rgSet_QC2)$Sample_Well %in% p1_failed_expQC)) |
                                                       (colData(stanford_meth_rgSet_QC2)$Sample_Plate == "Plt2" &
                                                          (colData(stanford_meth_rgSet_QC2)$Sample_Well %in% p2_failed_expQC)) |
                                                      (colData(stanford_meth_rgSet_QC2)$Sample_Plate == "Plt3" &
                                                      (colData(stanford_meth_rgSet_QC2)$Sample_Well %in% p3_failed_expQC)) |
                                                     (colData(stanford_meth_rgSet_QC2)$Sample_Plate == "Plt4" &
                                                       (colData(stanford_meth_rgSet_QC2)$Sample_Well %in% p4_failed_expQC))
                                                   ))
                                                   ]

stanford_meth_rgSet_QC3
# class: RGChannelSet 
# dim: 1051815 263 
# metadata(0):
#   assays(2): Green Red
# rownames(1051815): 1600101 1600111 ... 99810990 99810992
# rowData names(0):
#   colnames(263): Plt-1-C1-GIST Plt-1-D1-LMS ... Plt4-G2-GIST Plt4-H2-LM_old
# colData names(11): Sample_Name Sample_Well ... metadata ArrayTypes
# Annotation
# array: IlluminaHumanMethylationEPIC
# annotation: ilm10b4.hg19


# Save the dataset with the outlier samples removed
qsave(stanford_meth_rgSet_QC3, file=file.path(results_dir, "2_stanford_meth_rgSet_QC3.qs"), nthread=nthread)

## Calculate p values for all probes for the remaining samples
stanford_meth_pvals <- detectionP(stanford_meth_rgSet_QC3)

class(stanford_meth_pvals)
# [1] "matrix"

dim(stanford_meth_pvals)
# [1] 866091    263

# Identify genomic positions with p value > 0.01 ("failed" genomic positions)
failed <- stanford_meth_pvals > 0.01

# Calculate fraction of failed positions per sample
fractions_of_failed_positions <- colMeans(failed)

# Calculate a sum of failed positions per sample
sums_of_failed_positions <- colSums(failed)

# alternative way to calculate this:
# apply(stanford_meth_pvals, 2, sum)/nrow(stanford_meth_pvals)

# Save a file with the fraction of failed positions per sample
write.table(fractions_of_failed_positions , 
            file=file.path(qc_dir, "pval01_fractions_of_failed_positions_per_sample_263samples.tsv"), 
                        sep="\t ", na="NA", row.names=TRUE, col.names=TRUE)

# Save a file with the fraction of failed positions per sample
write.table(sums_of_failed_positions , 
            file=file.path(qc_dir, "pval01_sums_of_failed_positions_per_sample_263samples.tsv"), 
            sep="\t ", na="NA", row.names=TRUE, col.names=TRUE)

# How many genomic positions failed p value QC in >50% / >10% / >5% / >1% of samples? At least 1 sample (1/263 = 0.00380228)?
sum(rowMeans(failed) > 0.5) 
# [1] 602
sum(rowMeans(failed) > 0.1) 
# [1] 12240
sum(rowMeans(failed) > 0.05) 
# [1] 30482
sum(rowMeans(failed) > 0.01) 
# [1] 136804
sum(rowMeans(failed) > 0.003381) 
# [1] 269537


####
# Pre-Processing
####

print("Staring preprocessing...")

print("...Functional normalization")

## Normalization
## FIXME:: Is this supposed to be QC3?
stanford_methylSet <- preprocessFunnorm(stanford_meth_rgSet_QC3, ratioConvert = FALSE)

print("...minfiQC")

## Further QC
# Needed a methylSet to run this function
stanford_methylSet_QC <- minfiQC(stanford_methylSet, verbose = TRUE)

# Extract processed methylSet from returned list
stanford_methylSet <- stanford_methylSet_QC[[1]]

# Extract QC results dataframe from returned list
stanford_methylSet_QC <- stanford_methylSet_QC[[2]]; gc()

print("...Map to genome")

# Map to genome to create a GenomicMethylSet
stanford_gmSet <- mapToGenome(stanford_methylSet)

qsave(stanford_methylSet, file=file.path(results_dir, "3_stanford_methylSet.qs"), nthread=nthread)

qsave(stanford_gmSet, file=file.path(results_dir, "4_stanford_gmSet.qs"), nthread=nthread)

## Visualise what the QC3 data looks like before and after normalisation #####
pdf(file.path(qc_dir, "stanford_meth_rgSet_QC3.pdf"))
par(mfrow=c(1,2))

densityPlot(stanford_meth_rgSet_QC3, 
            sampGroups=colData(stanford_meth_rgSet_QC3)$Sample_Group, 
            main="Raw", 
            legend=FALSE, 
            pal=brewer.pal(8, "Set1"))
legend("topleft", 
      legend=levels(factor(colData(stanford_meth_rgSet_QC3)$Sample_Group)), 
      text.col=brewer.pal(8,"Set1"))

densityPlot(getBeta(stanford_methylSet), 
            sampGroups=colData(stanford_meth_rgSet_QC3)$Sample_Group,
            main="Normalized", 
            legend=FALSE, 
            pal = brewer.pal(8, "Set1"))
legend("topleft", 
       legend = levels(factor(colData(stanford_meth_rgSet_QC3)$Sample_Group)), 
       text.col=brewer.pal(8,"Set1"))

dev.off()

## Visualise what the QC2 data looks like before and after normalisation #####
pdf(file.path(qc_dir, "stanford_meth_rgSet_QC2"))
par(mfrow=c(1,2))
densityPlot(stanford_meth_rgSet_QC2, sampGroups=colData(stanford_meth_rgSet_QC2)$Sample_Group,main="Raw", legend=FALSE, pal = brewer.pal(8, "Set1"))
legend("topleft", legend = levels(factor(colData(stanford_meth_rgSet_QC2)$Sample_Group)), 
       text.col=brewer.pal(9,"Set1"))
densityPlot(getBeta(stanford_methylSet), sampGroups=colData(stanford_meth_rgSet_QC2)$Sample_Group,
            main="Normalized", legend=FALSE, pal = brewer.pal(8, "Set1"))
legend("topleft", legend = levels(factor(colData(stanford_meth_rgSet_QC2)$Sample_Group)), 
       text.col=brewer.pal(9,"Set1"))
dev.off()

rm(stanford_methylSet); gc()

print("...Convert to GenomicRatioSet")

## Convert to GenomicRatioSet
stanford_grSet <- ratioConvert(stanford_gmSet, what = "both")
rm(stanford_gmSet); gc()


## Check if sex chromosomes need to be removed
#### TODO:: determine if we want to exclude sex?
#### JP: I think yes. In the minfi 2014 paper they wrote: 
#"A multidimensional scaling plot of the methylation values reveals (not surprisingly) that sex is the biggest source of variability"

#### TODO:: add condition checking if there is more than one sex in the data
#table(getSex(stanford_LMS_grSet)$predictedSex) # 32 F 45 M

print("...annotate grSet")

# For now we will exclude sex as per Bioconductor Methylation analysis recommendation
# - Should better capture variability due to methylation status and not sex
stanford_grSet_annotations <- getAnnotation(stanford_grSet)
keep <- !(featureNames(stanford_grSet) %in%
  stanford_grSet_annotations$Name[
    stanford_grSet_annotations$chr %in% c("chrX", "chrY")]
  )
table(keep)
# keep
# FALSE   TRUE 
# 19627 846232 

stanford_grSet <- stanford_grSet[keep, ]

stanford_grSet
# class: GenomicRatioSet 
# dim: 846232 263 
# metadata(0):
#   assays(3): Beta M CN
# rownames(846232): cg14817997 cg26928153 ... cg07660283 cg09226288
# rowData names(0):
#   colnames(263): Plt-1-C1-GIST Plt-1-D1-LMS ... Plt4-G2-GIST Plt4-H2-LM_old
# colData names(16): Sample_Name Sample_Well ... mMed uMed
# Annotation
# array: IlluminaHumanMethylationEPIC
# annotation: ilm10b4.hg19
# Preprocessing
# Method: NA
# minfi version: NA
# Manifest version: NA

print("...filter poor quality probes")

## Filtering for Poor Quality Probes/Genomic positions ----

## Remove genomic positions that failed p value QC in at least 1 sample (1/263 = 0.00380228)
## Calculate p values for all probes for the remaining samples

# ensure probes are in the same order in the stanford_grSet and stanford_meth_pvals objects
stanford_meth_pvals <- stanford_meth_pvals[match(featureNames(stanford_grSet),
                                           rownames(stanford_meth_pvals)),] 

# remove any probes that have failed in one or more samples
keep <- rowSums(stanford_meth_pvals < 0.01) == ncol(stanford_grSet) 
table(keep)
stanford_grSetFlt <- stanford_grSet[keep,]
stanford_grSetFlt

# class: GenomicRatioSet 
# dim: 584901 263 
# metadata(0):
#   assays(3): Beta M CN
# rownames(584901): cg26928153 cg16269199 ... cg07660283 cg09226288
# rowData names(0):
#   colnames(263): Plt-1-C1-GIST Plt-1-D1-LMS ... Plt4-G2-GIST Plt4-H2-LM_old
# colData names(16): Sample_Name Sample_Well ... mMed uMed
# Annotation
# array: IlluminaHumanMethylationEPIC
# annotation: ilm10b4.hg19
# Preprocessing
# Method: NA
# minfi version: NA
# Manifest version: NA

rm(stanford_grSet); gc()

## SNP Correction ----

print("...SNP correction")

# Get SNP info dataframe
stanford_SNPs <- getSnpInfo(stanford_grSetFlt)
length(grep(stanford_SNPs$CpG_rs, pattern = "rs*")) # Number of CpG SNPs
# [1] 17370
length(grep(stanford_SNPs$SBE_rs, pattern = "rs*")) # Number of SBE SNPs
# [1] 9071

# Add SNP info to grSet object
stanford_grSetFlt <- addSnpInfo(stanford_grSetFlt)

# Drop probes containing a SNP at CpG interrogation or at single NT extension
stanford_grSetFlt <- dropLociWithSnps(stanford_grSetFlt, snps=c("SBE","CpG"), maf=0)

stanford_grSetFlt
# class: GenomicRatioSet 
# dim: 567070 263 
# metadata(0):
#   assays(3): Beta M CN
# rownames(567070): cg26928153 cg16269199 ... cg07660283 cg09226288
# rowData names(6): Probe_rs Probe_maf ... SBE_rs SBE_maf
# colnames(263): Plt-1-C1-GIST Plt-1-D1-LMS ... Plt4-G2-GIST Plt4-H2-LM_old
# colData names(16): Sample_Name Sample_Well ... mMed uMed
# Annotation
# array: IlluminaHumanMethylationEPIC
# annotation: ilm10b4.hg19
# Preprocessing
# Method: NA
# minfi version: NA
# Manifest version: NA


## Cross-reactivity Correction ----

print("...Cross-reactive probe correction")

# Get list of EPIC xreactive probes
EPIC_xreactive_probes <- maxprobes::xreactive_probes()

# Filter out xreactive probes 
keep <- !(featureNames(stanford_grSetFlt) %in% EPIC_xreactive_probes)
table(keep)

stanford_grSetFlt <- stanford_grSetFlt[keep, ]
stanford_grSetFlt
# class: GenomicRatioSet 
# dim: 532745 263 
# metadata(0):
#   assays(3): Beta M CN
# rownames(532745): cg26928153 cg16269199 ... cg19565306 cg09226288
# rowData names(6): Probe_rs Probe_maf ... SBE_rs SBE_maf
# colnames(263): Plt-1-C1-GIST Plt-1-D1-LMS ... Plt4-G2-GIST Plt4-H2-LM_old
# colData names(16): Sample_Name Sample_Well ... mMed uMed
# Annotation
# array: IlluminaHumanMethylationEPIC
# annotation: ilm10b4.hg19
# Preprocessing
# Method: NA
# minfi version: NA
# Manifest version: NA

qsave(stanford_grSetFlt, 
      file=file.path(results_dir, "6_stanford_grSet_filtered_preprocessed_263samples.qs"), 
      nthread=nthread)

####
# Subset data
####

print('Staring data subsetting...')

## Filtering to samples of interest

# All LMS samples
print("...LMS samples")

stanford_LMS_grSetFlt <- stanford_grSetFlt[, (colData(stanford_grSetFlt)$Sample_Group == "LMS" |
                                          colData(stanford_grSetFlt)$Sample_Group == "uLMS")]

qsave(stanford_LMS_grSetFlt, file = file.path(results_dir, "6a_stanford_LMS_grSetFlt.qs"), nthread=nthread)



# LM samples

print("...LM samples")

## FIXME:: Is the space in LM old a typo?
stanford_LM_grSetFlt <- stanford_grSetFlt[, (colData(stanford_grSetFlt)$Sample_Group == "LM old " |
                                         colData(stanford_grSetFlt)$Sample_Group == "LM new")]

qsave(stanford_LM_grSetFlt, file = file.path(results_dir, "6b_stanford_LM_grSetFlt.qs"), nthread=nthread)
rm(stanford_LM_gmSet); gc()

# GIST samples
print("...GIST samples")

stanford_GIST_grSetFlt <- stanford_grSetFlt[, (colData(stanford_grSetFlt)$Sample_Group == "GIST" |
                                           colData(stanford_grSetFlt)$Sample_Group == "miniGIST")]

qsave(stanford_GIST_grSetFlt, file = file.path(results_dir, "6c_stanford_GIST_grSetFlt.qs"), nthread=nthread)


# All DTF samples
print("...DTF samples")

stanford_DTF_grSetFlt <- stanford_grSetFlt[, (colData(stanford_grSetFlt)$Sample_Group == "DTF" |
                                          colData(stanford_grSetFlt)$Sample_Group == "DTF_het")]

qsave(stanford_DTF_grSetFlt, file = file.path(results_dir, "6d_stanford_DTF_grSetFlt.qs"), nthread=nthread )

# DTF_heterogeneity_samples
print("...DTF heterogeneity samples")
stanford_DTFhet_grSetFlt <- stanford_grSetFlt[, colData(stanford_grSetFlt)$Sample_Group == "DTF_het"]

qsave(stanford_DTFhet_grSetFlt, file = file.path(results_dir, "6e_stanford_DTFhet_grSetFlt.qs"), nthread=nthread)



####
# Differential Methylation Analysis of Regions
# - Using M value due to better distributional properties
# - Don't have a control group to compare against; can we get historic FFPE for health tissue?
####

print("Starting differential methylation analysis...")

## Currently have no control to run DMR against, going to use feature/sample Beta values

# Extract M values (feature x sample)
stanford_263samples_M_values <- getM(stanford_grSetFlt)

# Write to .csv file for SNF
write.csv(stanford_263samples_M_values, file = file.path(results_dir, "7a_stanford_263samples_feature_M_matrix.csv"))

print("...Extracting beta values")

# Extract beta values
stanford_263samples_bVals <- getBeta(stanford_grSetFlt)

stanford_LMS_bVals <- getBeta(stanford_LMS_grSetFlt)

stanford_GIST_bVals <- getBeta(stanford_GIST_grSetFlt)

stanford_DTF_bVals <- getBeta(stanford_DTF_grSetFlt)

stanford_DTFhet_bVals <- getBeta(stanford_DTFhet_grSetFlt)

print("...Writing to csv")

# Write to .csv file for SNF
write.csv(stanford_263samples_bVals, file = file.path(results_dir, "7b_stanford_263samples_feature_beta_matrix.csv"))
write.csv(stanford_LMS_bVals, file = file.path(results_dir, "7d_stanford_LMS_feature_beta_matrix.csv"))
write.csv(stanford_GIST_bVals, file = file.path(results_dir, "7e_stanford_GIST_feature_beta_matrix.csv"))
write.csv(stanford_DTF_bVals, file = file.path(results_dir, "7f_stanford_DTF_feature_beta_matrix.csv"))
write.csv(stanford_DTFhet_bVals, file = file.path(results_dir, "7g_stanford_DTFhet_feature_beta_matrix.csv"))

print("...Filter to top 1k most variant features")

## Filter top 1000 genomic positions with the highest variance of beta values to use for hierarchical clustering and PCA
vars <- rowVars(stanford_263samples_bVals)
d <- stanford_263samples_bVals[rev(order(vars))[seq_len(1000)], ]
dim(d)

# Save matrix with 1000 top variable CpGs
write.csv(d, file=file.path(results_dir, "7c_stanford_263samples_feature_top1000_beta_matrix.csv"))

## Filter top 1000 genomic positions with the highest variance of beta values to use for hierarchical clustering and PCA for LMS samples
vars <- rowVars(stanford_LMS_bVals)
d <- stanford_LMS_bVals[rev(order(vars))[seq_len(1000)], ]
dim(d)
gc()

# Save matrix with 1000 top variable CpGs
write.csv(d, file=file.path(results_dir, "7h_stanford_LMS_feature_top1000_beta_matrix.csv"))

## Filter top 1000 genomic positions with the highest variance of beta values to use for hierarchical clustering and PCA for GIST samples
vars <- rowVars(stanford_GIST_bVals)
d <- stanford_GIST_bVals[rev(order(vars))[1:1000], ]
dim(d)
gc()

# Save matrix with 1000 top variable CpGs
write.csv(d, file=file.path(results_dir, "7i_stanford_GIST_feature_top1000_beta_matrix.csv"))

## Filter top 1000 genomic positions with the highest variance of beta values to use for hierarchical clustering and PCA for DTF samples
vars <- rowVars(stanford_DTF_bVals)
d <- stanford_DTF_bVals[rev(order(vars))[1:1000], ]
dim(d)
gc()

# Save matrix with 1000 top variable CpGs
write.csv(d, file=file.path(results_dir, "7j_stanford_DTF_feature_top1000_beta_matrix.csv"))

## Filter top 1000 genomic positions with the highest variance of beta values to use for hierarchical clustering and PCA for DTF_het samples
vars <- rowVars(stanford_DTFhet_bVals)
d <- stanford_DTFhet_bVals[rev(order(vars))[1:1000], ]
dim(d)
gc()

# Save matrix with 1000 top variable CpGs
write.csv(d, file = file.path(results_dir, "7k_stanford_DTFhet_feature_top1000_beta_matrix.csv"))

rm(d); gc()

print("...Collapse adjacent CpG sites")

# Collapse adjacent CpG sites
stanford_grSetFlt_M_collapsed <- cpgCollapse(stanford_grSetFlt, what = c('M'))
qsave(stanford_grSetFlt_M_collapsed, file = file.path(results_dir, '8a_stanford_263samples_grSet_M_collapsed.qs'), nthread=nthread)

stanford_grSetFlt_Beta_collapsed <- cpgCollapse(stanford_grSetFlt, what = c('Beta'))
qsave(stanford_grSetFlt_Beta_collapsed, file = file.path(results_dir, '8b_stanford_263samples_grSet_Beta_collapsed.qs'), nthread=nthread)


# Mapping data from collapsed sites
block_to_cpg_mapping <- stanford_grSetFlt_M_collapsed[[2]]$indexes
block_to_cpg_mapping <- stanford_grSetFlt_Beta_collapsed[[2]]$indexes

# Extract feature names
features_M <- rownames(stanford_263samples_M_values)
features_Beta <- rownames(stanford_263samples_bVals)


# Get the maximum number of CpG in any block
max_cpgs_per_block <- max(vapply(block_to_cpg_mapping, length, FUN.VALUE=numeric(1)))

# Rename the cpg indexes to cpg names and fill NAs up to max number of cpg
block_to_cpg_mapping_M <- lapply(block_to_cpg_mapping, 
                                 function(block) 
                                    c(features_M[block], rep(NA, max_cpgs_per_block - length(block))))
block_to_cpg_mapping_Beta <- lapply(block_to_cpg_mapping, 
                                    function(block) 
                                        c(features_Beta[block], rep(NA, max_cpgs_per_block - length(block))))

print("...Merging collapsed data and writing mapping files")

# Bind the list into a data.frame
block_to_cpg_mapping_M_df <- do.call(cbind, block_to_cpg_mapping_M)
block_to_cpg_mapping_Beta_df <- do.call(cbind, block_to_cpg_mapping_Beta)

# Write mapping to .csv file
write.csv(block_to_cpg_mapping_M_df, file= paste0(results_dir, '/9a_block_to_cpg_mapping_M.csv'))
write.csv(block_to_cpg_mapping_Beta_df, file= paste0(results_dir, '/9a_block_to_cpg_mapping_Beta.csv'))


print("...Writing clustering M value data to csv")

# Write to .csv file for SNF
write.csv(getM(stanford_grSetFlt_M_collapsed[[1]]), file = paste0(results_dir, "/9b_stanford_LMS_sample_v_block_M.csv"))

## Filter top 1000 blocks with the highest variance of beta values to use for hierarchical clustering and PCA
vars <- rowVars(block_to_cpg_mapping_Beta_df)
d6 <- block_to_cpg_mapping_Beta_df[rev(order(vars))[1:1000], ]
dim(d6)

# Save matrix with 1000 top variable CpGs
write.csv(d6, file=paste0(results_dir, "/7c_stanford_263samples_feature_top1000_beta_matrix.csv"))

print("Finishing script...")
print("...sessionInfo:")
print(sessionInfo())

###
# END
###