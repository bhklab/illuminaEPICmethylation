###
# Illumina Epic LMS Methylation Analysis with minfi
# by Christopher Eeles
#
# Adapted from:
# "A cross-package Bioconductor workflow for analysing methylation array data" from Bioconductor
# - Available at: https://urldefense.com/v3/__https://www.bioconductor.org/packages/devel/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html__;!!CjcC7IQ!eq7O-K4MIWeffhyGDxrjE9ueXefAlaXS84G55CewG_nAkXP3t1lHw35OToXkg9BOz-h4EfrbL_I$ 
# "BioC2014 - Minfi tutorial"
# - Available at: https://urldefense.com/v3/__https://www.bioconductor.org/help/course-materials/2014/BioC2014/minfi_BioC2014.pdf__;!!CjcC7IQ!eq7O-K4MIWeffhyGDxrjE9ueXefAlaXS84G55CewG_nAkXP3t1lHw35OToXkg9BOz-h48Vi55BI$ 
# ""
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

# Clear any objects in environment to save on RAM
rm(list=ls())

####
# Data Import
####

# Configure path to data
data_dir <- '~/Desktop/Methylation_analysis_April2020/EPIC_Methylation_data_2019'
metadata_dir <- '~/Desktop/Methylation_analysis_April2020/labels_copy'
results_dir <- '~/Desktop/Methylation_analysis_April2020/results'
qc_dir <- '~/Desktop/Methylation_analysis_April2020/qc_results'

## Read data by plate


# Configure plate paths
p1_dir <- paste0(data_dir, "/Plate1")
p2_dir <- paste0(data_dir, "/Plate2")
p3_dir <- paste0(data_dir, "/Plate3")

# Create methylation array sheets
p1_sheet <- read.metharray.sheet(p1_dir)
p2_sheet <- read.metharray.sheet(p2_dir)
p3_sheet <- read.metharray.sheet(p3_dir)

# Import plate labels
p1_labels <- read.csv(paste0(metadata_dir, "/plate1_labels.csv"), stringsAsFactors = FALSE)
p2_labels <- read.csv(paste0(metadata_dir, "/plate2_labels.csv"), stringsAsFactors = FALSE)
p3_labels <- read.csv(paste0(metadata_dir, "/plate3_labels.csv"), stringsAsFactors = FALSE)

# Label sample groups
p1_sheet$Sample_Group <- p1_labels$Group
p2_sheet$Sample_Group <- p2_labels$Group
p3_sheet$Sample_Group <- p3_labels$Group

# Read data into Extended rgSets
p1_RG <- read.metharray.exp(targets = p1_sheet, extended = TRUE)
p2_RG <- read.metharray.exp(targets = p2_sheet, extended = TRUE)
p3_RG <- read.metharray.exp(targets = p3_sheet, extended = TRUE)

# Read data into rgSets
# p1_RG <- read.metharray.exp(targets = p1_sheet)
# p2_RG <- read.metharray.exp(targets = p2_sheet)
# p3_RG <- read.metharray.exp(targets = p3_sheet)

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

## Merge the plates into one rgSet
p12_RG <- combineArrays(p1_RG, p2_RG, outType = "IlluminaHumanMethylationEPIC")
stanford_meth_rgSet <- combineArrays(p12_RG, p3_RG, outType = "IlluminaHumanMethylationEPIC")

## JP: check content of the extended stanford_meth_rgSet:
stanford_meth_rgSet
# class: RGChannelSetExtended 
# dim: 1051815 304 
# metadata(0):
#   assays(5): Green Red GreenSD RedSD NBeads
# rownames(1051815): 1600101 1600111 ... 99810990 99810992
# rowData names(0):
#   colnames(304): Plt-1-A1-Control Plt-1-B1-DTF_het ... Plt4-G2-GIST Plt4-H2-LM
# colData names(11): Sample_Name Sample_Well ... metadata ArrayTypes
# Annotation
# array: IlluminaHumanMethylationEPIC
# annotation: ilm10b4.hg19

## JP: check content of the stanford_meth_rgSet:
# stanford_meth_rgSet
# class: RGChannelSet 
# dim: 1051815 304 
# metadata(0):
#   assays(2): Green Red
# rownames(1051815): 1600101 1600111 ... 99810990 99810992
# rowData names(0):
#   colnames(304): Plt-1-A1-Control Plt-1-B1-DTF_het ... Plt4-G2-GIST Plt4-H2-LM
# colData names(11): Sample_Name Sample_Well ... metadata ArrayTypes
# Annotation
# array: IlluminaHumanMethylationEPIC
# annotation: ilm10b4.hg19


## Save the extended stanford_meth_rgSet object
save(stanford_meth_rgSet, file = paste0(results_dir, "/1_stanford_meth_rgSetExtended_raw.RData"))
rm(p1_RG, p2_RG, p3_RG, p12_RG)
gc()

## Save stanford_meth_rgSet object
# save(stanford_meth_rgSet, file = paste0(results_dir, "/1_stanford_meth_rgSet_raw.RData"))
# rm(p1_RG, p2_RG, p3_RG, p12_RG)
# gc()


####
# Quality Control
####

## Remove probes failing experimental QC

# Per plate failed sample wells (16 wells)
# Based on the QC plots from the facility - samples with < 800k CpGs with p value < 0.05
# p1_failed_expQC <- c("E1", "D2", "A3", "G3", "H3", "A4", "B5", "H5", "A7", "B12")
# p2_failed_expQC <- c("C4", "D5", "G5", "D6", "F6", "G6")
# p3_failed_expQC <- c()
# p4_failed_expQC <- c()

# Subset to rgSet not containing the failed sample wells
#stanford_meth_rgSet <- stanford_meth_rgSet[,
#  which(!(
#    (colData(stanford_meth_rgSet)$Sample_Plate == "Plt-1" &
#      (colData(stanford_meth_rgSet)$Sample_Well %in% p1_failed_expQC)) |
#    (colData(stanford_meth_rgSet)$Sample_Plate == "Plt2" &
#      (colData(stanford_meth_rgSet)$Sample_Well %in% p2_failed_expQC)) # |
    # (colData(stanford_meth_rgSet)$Sample_Plate == "Plt3" &
    #   (colData(stanford_meth_rgSet)$Sample_Well %in% p3_failed_expQC)) |
    # (colData(stanford_meth_rgSet)$Sample_Plate == "Plt4" &
    #   (colData(stanford_meth_rgSet)$Sample_Well %in% p4_failed_expQC))
#  ))
# ]

## JP: check content of the stanford_meth_rgSet after subsetting:
#stanford_meth_rgSet
# class: RGChannelSet 
# dim: 1051815 288 
# metadata(0):
#   assays(2): Green Red
# rownames(1051815): 1600101 1600111 ... 99810990 99810992
# rowData names(0):
#   colnames(288): Plt-1-A1-Control Plt-1-B1-DTF_het ... Plt4-G2-GIST Plt4-H2-LM_old
# colData names(11): Sample_Name Sample_Well ... metadata ArrayTypes
# Annotation
# array: IlluminaHumanMethylationEPIC
# annotation: ilm10b4.hg19

## Calculate pvals using detectionP() function

# Default usage: detectionP(rgSet, type = "m+u")
# The m+u method compares the total DNA signal (Methylated + Unmethylated) 
# for each position to the background signal level.

# A detection p-value is returned FOR EVERY GENOMIC POSITION in every sample. 
# (that's why it returns a matrix with ~866k rows, and not ~1mln rows) 
# Small p-values indicate a good position. 
# Positions with non-significant p-values (typically >0.01) should not be trusted.

# Most of the pvals are 0, but this is also the case in the Bioconductor guide
# stanford_meth_pvals <- detectionP(stanford_meth_rgSet)

# class(stanford_meth_pvals)
# [1] "matrix"

# dim(stanford_meth_pvals)
# [1] 866091    288

## Chris: Remove probes failing pval QC - JP comment: this formula calculates mean p value per sample and removes samples with mean p value > 0.05 - not sure if this is a valid approach
# stanford_meth_rgSet_QC1 <- stanford_meth_rgSet[, colMeans(stanford_meth_pvals) < 0.05 ]

## Dimenstions of stanford_meth_rgSet_QC1 when < 0.05
# dim(stanford_meth_rgSet_QC1)
# [1] 1051815     286

# rm(stanford_meth_rgSet) # This will save lots of RAM
# gc()
# save(stanford_meth_rgSet_QC1, file = paste0(results_dir, "/2_stanford_meth_rgSet_QC1.RData"))

## Filter out samples removed in previous step
# stanford_meth_pvals <- stanford_meth_pvals[, colMeans(stanford_meth_pvals) < 0.01]

## Generate QC Report
# qcReport(
#   stanford_meth_rgSet_QC1,
#   sampNames = sampleNames(stanford_meth_rgSet_QC1),
#   sampGroups = colData(stanford_meth_rgSet_QC1)$Sample_Group,
#   pdf = paste0(qc_dir, "/stanford_meth_QC1_report.pdf")
# )

## Filter further based on QC Report
#### TO-DO:: Determine how to interpret beta plots
# - Can't find good documentation on interpreting QC plots

# # Sample wells that failed visual QC in report
# p1_failed_QC2 <- c()
# p2_failed_QC2 <- c()
# p3_failed_QC2 <- c()
# p4_failed_QC2 <- c()

# # Subset to all samples not failing QC2
# stanford_meth_rgSet_QC2 <- stanford_meth_rgSet_QC1[,
#   which(!(
#     (colData(stanford_meth_rgSet)$Sample_Plate == "Plt-1" &
#       (colData(stanford_meth_rgSet)$Sample_Well %in% p1_failed_QC2)) |
#     (colData(stanford_meth_rgSet)$Sample_Plate == "Plt2" &
#       (colData(stanford_meth_rgSet)$Sample_Well %in% p2_failed_QC2)) # |
#     # (colData(stanford_meth_rgSet)$Sample_Plate == "Plt3" &
#     #   (colData(stanford_meth_rgSet)$Sample_Well %in% p3_failed_QC2)) |
#     # (colData(stanford_meth_rgSet)$Sample_Plate == "Plt4" &
#     #   (colData(stanford_meth_rgSet)$Sample_Well %in% p4_failed_QC2))
#   ))
# ]
# rm(stanford_meth_rgSet_QC1)
# save(stanford_meth_rgSet_QC2, file = "stanford_meth_rgSet_QC2.RData")


## Remove this when QC2 is completed
# stanford_meth_rgSet_QC2 <- stanford_meth_rgSet_QC1
# rm(stanford_meth_rgSet_QC1)
# gc()


#####################################################################################################################################################################################
#####################################################################################################################################################################################
#### JP: alternative QC -exclude outlier samples and calculate p values #############################################################################################################
#####################################################################################################################################################################################
#####################################################################################################################################################################################


#### Illumina's recommendation for the QC of FFPE EPIC arrays: 
#### keep samples with >90% detected CpGs at p < 0.05

## Calculate p values for all probes for all samples
stanford_meth_pvals <- detectionP(stanford_meth_rgSet)

class(stanford_meth_pvals)
# [1] "matrix"

dim(stanford_meth_pvals)
# [1] 866091    304


# Identify genomic positions with p value < 0.05 (detected genomic positions)
detected <- stanford_meth_pvals < 0.05

# Calculate fraction of failed positions per sample
fractions_of_detected_positions <- colMeans(detected)

# Calculate a sum of failed positions per sample
sums_of_detected_positions <- colSums(detected)

# alternative way to calculate this:
# apply(stanford_meth_pvals, 2, sum)/nrow(stanford_meth_pvals)

# Save a file with the fraction of failed positions per sample
write.table(fractions_of_detected_positions , file = paste0(qc_dir, "/fractions_of_detected_positions_p0_05_per_sample_304samples.txt"),sep = "\t ", na = "NA", row.names = TRUE, col.names = TRUE)

# Save a file with the number of detected positions per sample
write.table(sums_of_detected_positions , file = paste0(qc_dir, "/sums_of_detected_positions_p0_05_per_sample_304samples.txt"),sep = "\t ", na = "NA", row.names = TRUE, col.names = TRUE)

# How many samples have > 90% of CpGs detected at  p value < 0.05?
sum(colMeans(detected) > 0.9) 
# [1] 294

# Which samples have < 90% CpGs detected at p value < 0.05 and should be excluded?
which(colMeans(detected) < 0.9)
# Plt-1-D2-uLMS Plt-1-B5-GIST  Plt-1-A7-LMS   Plt2-D5-LMS   Plt2-F5-LMS   Plt2-D6-LMS  Plt2-F6-uLMS  Plt2-G6-uLMS  Plt3-E6-uLMS  Plt3-C7-uLMS 
# 12            34            49           132           134           140           142           143           237           243 


## Remove 10 samples that have < 90% CpGs detected at p value < 0.05:

p1_failed_expQC <- c("A7", "B5", "D2")
p2_failed_expQC <- c("D5", "D6", "F5", "F6", "G6")
p3_failed_expQC <- c("C7", "E6")
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
# class: RGChannelSetExtended
# dim: 1051815 294 
# metadata(0):
#   assays(5): Green Red GreenSD RedSD NBeads
# rownames(1051815): 1600101 1600111 ... 99810990 99810992
# rowData names(0):
#   colnames(294): Plt-1-A1-Control Plt-1-B1-DTF_het ... Plt4-G2-GIST Plt4-H2-LM
# colData names(11): Sample_Name Sample_Well ... metadata ArrayTypes
# Annotation
# array: IlluminaHumanMethylationEPIC
# annotation: ilm10b4.hg19

# stanford_meth_rgSet_QC1
# class: RGChannelSet 
# dim: 1051815 294 
# metadata(0):
#   assays(2): Green Red
# rownames(1051815): 1600101 1600111 ... 99810990 99810992
# rowData names(0):
#   colnames(294): Plt-1-A1-Control Plt-1-B1-DTF_het ... Plt4-G2-GIST Plt4-H2-LM
# colData names(11): Sample_Name Sample_Well ... metadata ArrayTypes
# Annotation
# array: IlluminaHumanMethylationEPIC
# annotation: ilm10b4.hg19


# Save the extended dataset with the outlier samples removed
save(stanford_meth_rgSet_QC1, file = paste0(results_dir, "/2_stanford_meth_rgSetExtended_QC1.RData"))

# Save the dataset with the outlier samples removed
# save(stanford_meth_rgSet_QC1, file = paste0(results_dir, "/2_stanford_meth_rgSet_QC1.RData"))

#### JP: Check bisulfite conversion rate using bscon function from the waterMelon package

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("wateRmelon")

library(wateRmelon)

bs <- bscon(stanford_meth_rgSet_QC1)
head(bs)
range(bs)
# [1] 57.70437 97.17750
median(bs)
# [1] 92.62639
mean(bs)
# [1] 92.11485

# Save a file with bisulfite conversion rate for each sample
write.csv(bs, file=paste0(qc_dir, "/Bisulfite_conversion_rate_QC1.csv"))

# How many samples have < 90% bisulfite conversion efficiency?
sum(bs < 90)
# [1] 28

# Which samples have < 90% bisulfite conversion efficiency?
which(bs < 90)
# Plt-1-D1-LMS     Plt-1-C2-GIST     Plt-1-E2-ULMS     Plt-1-H2-GIST      Plt-1-H3-LMS      Plt-1-A4-LMS      Plt-1-H5-LMS      Plt-1-A6-LMS 
# 4                11                12                15                23                24                38                39 
# Plt-1-G7-LMS    Plt-1-E10-GIST    Plt-1-H10-GIST    Plt-1-B11-ULMS Plt-1-G11-Control       Plt2-C2-LMS       Plt2-D2-LMS       Plt2-D3-LMS 
# 52                74                77                79                84               104               105               113 
# Plt2-H4-LM        Plt2-H5-LM        Plt2-B8-LM      Plt2-E10-LMS      Plt3-D2-ULMS        Plt3-E2-LM   Plt3-H3-DTF_het       Plt3-A4-LMS 
# 125               131               146               165               196               197               208               209 
# Plt3-C8-ULMS     Plt3-B10-ULMS     Plt3-A12-ULMS   Plt4-B1-DTF_het 
# 241               256               271               280 


# # How many samples have < 85% bisulfite conversion efficiency?
# sum(bs < 85)
# # [1] 11
# 
# # Which samples have < 85% bisulfite conversion efficiency?
# which(bs < 85)
# # Plt-1-C2-GIST    Plt-1-H3-LMS    Plt-1-H5-LMS    Plt-1-G7-LMS  Plt-1-E10-GIST     Plt2-D2-LMS      Plt2-H4-LM      Plt2-H5-LM    Plt3-D2-ULMS 
# # 11              23              38              52              74             105             125             131             196 
# # Plt3-H3-DTF_het    Plt3-C8-ULMS 
# # 208             241 
# 
# # How many samples have < 80% bisulfite conversion efficiency?
# sum(bs < 80)
# # [1] 4
# 
# # Which samples have < 80% bisulfite conversion efficiency?
# which(bs < 80)
# # Plt2-H4-LM    Plt3-D2-ULMS Plt3-H3-DTF_het    Plt3-C8-ULMS 
# # 125             196             208             241 

## Remove 28 samples that have < 90% bisulfite conversion efficiency:

p1_failed_expQC <- c("D1", "C2", "E2", "H2", "H3", "A4", "H5", "A6", "G7", "E10", "H10", "B11", "G11")
p2_failed_expQC <- c("C2", "D2", "D3", "H4", "H5", "B8", "E10")
p3_failed_expQC <- c("D2", "E2", "H3", "A4", "C8", "B10", "A12")
p4_failed_expQC <- c("B1")

# ## Remove 11 samples that have < 85% bisulfite conversion efficiency:
# 
# p1_failed_expQC <- c("C2", "H3", "H5", "G7", "E10")
# p2_failed_expQC <- c("D2", "H4", "H5")
# p3_failed_expQC <- c("D2", "H3", "C8")
# p4_failed_expQC <- c()
# 
# ## Remove 4 samples that have < 80% bisulfite conversion efficiency:
# 
# p1_failed_expQC <- c()
# p2_failed_expQC <- c("H4")
# p3_failed_expQC <- c("D2", "H3", "C8")
# p4_failed_expQC <- c()

stanford_meth_rgSet_QC2 <- stanford_meth_rgSet_QC1[,
                                               which(!(
                                                 (colData(stanford_meth_rgSet_QC1)$Sample_Plate == "Plt-1" &
                                                    (colData(stanford_meth_rgSet_QC1)$Sample_Well %in% p1_failed_expQC)) |
                                                   (colData(stanford_meth_rgSet_QC1)$Sample_Plate == "Plt2" &
                                                      (colData(stanford_meth_rgSet_QC1)$Sample_Well %in% p2_failed_expQC)) |
                                                   (colData(stanford_meth_rgSet_QC1)$Sample_Plate == "Plt3" &
                                                      (colData(stanford_meth_rgSet_QC1)$Sample_Well %in% p3_failed_expQC)) |
                                                   (colData(stanford_meth_rgSet_QC1)$Sample_Plate == "Plt4" &
                                                     (colData(stanford_meth_rgSet_QC1)$Sample_Well %in% p4_failed_expQC))
                                               ))
                                               ]

stanford_meth_rgSet_QC2
# class: RGChannelSetExtended 
# dim: 1051815 266 
# metadata(0):
#   assays(5): Green Red GreenSD RedSD NBeads
# rownames(1051815): 1600101 1600111 ... 99810990 99810992
# rowData names(0):
#   colnames(266): Plt-1-A1-Control Plt-1-B1-DTF_het ... Plt4-G2-GIST Plt4-H2-LM
# colData names(11): Sample_Name Sample_Well ... metadata ArrayTypes
# Annotation
# array: IlluminaHumanMethylationEPIC
# annotation: ilm10b4.hg19

# stanford_meth_rgSet_QC2
# class: RGChannelSet 
# dim: 1051815 290 
# metadata(0):
#   assays(2): Green Red
# rownames(1051815): 1600101 1600111 ... 99810990 99810992
# rowData names(0):
#   colnames(290): Plt-1-A1-Control Plt-1-B1-DTF_het ... Plt4-G2-GIST Plt4-H2-LM
# colData names(11): Sample_Name Sample_Well ... metadata ArrayTypes
# Annotation
# array: IlluminaHumanMethylationEPIC
# annotation: ilm10b4.hg19

# Save the dataset with the outlier samples with poor bisulfite conversion efficiency removed
# save(stanford_meth_rgSet_QC2, file = paste0(results_dir, "/2_stanford_meth_rgSet_BC80_QC2.RData"))
# save(stanford_meth_rgSet_QC2, file = paste0(results_dir, "/2_stanford_meth_rgSetExtended_BC85_QC2.RData"))
save(stanford_meth_rgSet_QC2, file = paste0(results_dir, "/2_stanford_meth_rgSetExtended_BC90_QC2.RData"))

## Generate QC Report on samples from QC2
qcReport(
  stanford_meth_rgSet_QC2,
  sampNames = sampleNames(stanford_meth_rgSet_QC2),
  sampGroups = colData(stanford_meth_rgSet_QC2)$Sample_Group,
  pdf = paste0(qc_dir, "/stanford_meth_QC2_report.pdf")
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



## shinyMethyl is particularly useful to visualize all plots at the same time in an interactive fashion.
library(shinyMethyl)
myShinyMethylSet = shinySummarize(stanford_meth_rgSet_QC2)
runShinyMethyl(myShinyMethylSet)

## Based on the QC plot from shinyMethyl, remove 1 sample with log2 M or U below 9 (in the the dataset filtered by bisulfite conversion >80%)
p1_failed_expQC <- c("H3")
p2_failed_expQC <- c()
p3_failed_expQC <- c()
p4_failed_expQC <- c()

stanford_meth_rgSet_QC3 <- stanford_meth_rgSet_QC2[,
                                                   which(!(
                                                     (colData(stanford_meth_rgSet_QC2)$Sample_Plate == "Plt-1" &
                                                        (colData(stanford_meth_rgSet_QC2)$Sample_Well %in% p1_failed_expQC)) |
                                                       (colData(stanford_meth_rgSet_QC2)$Sample_Plate == "Plt2" &
                                                          (colData(stanford_meth_rgSet_QC2)$Sample_Well %in% p2_failed_expQC)) |
                                                       (colData(stanford_meth_rgSet_QC2)$Sample_Plate == "Plt3" &
                                                          (colData(stanford_meth_rgSet_QC2)$Sample_Well %in% p3_failed_expQC)) # |
                                                      # (colData(stanford_meth_rgSet_QC2)$Sample_Plate == "Plt4" &
                                                      #    (colData(stanford_meth_rgSet_QC2)$Sample_Well %in% p4_failed_expQC))
                                                   ))
                                                   ]

stanford_meth_rgSet_QC3
# class: RGChannelSet 
# dim: 1051815 289 
# metadata(0):
#   assays(2): Green Red
# rownames(1051815): 1600101 1600111 ... 99810990 99810992
# rowData names(0):
#   colnames(289): Plt-1-A1-Control Plt-1-B1-DTF_het ... Plt4-G2-GIST Plt4-H2-LM
# colData names(11): Sample_Name Sample_Well ... metadata ArrayTypes
# Annotation
# array: IlluminaHumanMethylationEPIC
# annotation: ilm10b4.hg19

# Save the dataset QC3 with the outlier samples removed
save(stanford_meth_rgSet_QC3, file = paste0(results_dir, "/2_stanford_meth_rgSet_logMU9_QC3.RData"))


## Run shinyMethyl on the filtered dataset
library(shinyMethyl)
myShinyMethylSet = shinySummarize(stanford_meth_rgSet_QC3)
runShinyMethyl(myShinyMethylSet)




###### REMOVE TECHNICAL REPLICATES OR SAMPLES FOR HETEROGENEITY STUDIES ######

## If no filtering based on shinyMethyl QC plot, then:
stanford_meth_rgSet_QC3 <- stanford_meth_rgSet_QC2

## Remove 16 Jurkat DNA samples (these are just controls for assay performance)
## I think fully methylated Jurkat DNA may affect normalization?
p1_failed_expQC <- c("B2", "D5", "F9", "A1", "C4", "E8", "G11")
p2_failed_expQC <- c("F3", "E6", "C8", "A12", "E4", "D7", "B11")
p3_failed_expQC <- c("B7")
p4_failed_expQC <- c("A1")

stanford_meth_rgSet_QC4 <- stanford_meth_rgSet_QC3[,
                                                   which(!(
                                                     (colData(stanford_meth_rgSet_QC3)$Sample_Plate == "Plt-1" &
                                                        (colData(stanford_meth_rgSet_QC3)$Sample_Well %in% p1_failed_expQC)) |
                                                       (colData(stanford_meth_rgSet_QC3)$Sample_Plate == "Plt2" &
                                                          (colData(stanford_meth_rgSet_QC3)$Sample_Well %in% p2_failed_expQC)) |
                                                       (colData(stanford_meth_rgSet_QC3)$Sample_Plate == "Plt3" &
                                                          (colData(stanford_meth_rgSet_QC3)$Sample_Well %in% p3_failed_expQC)) |
                                                       (colData(stanford_meth_rgSet_QC3)$Sample_Plate == "Plt4" &
                                                          (colData(stanford_meth_rgSet_QC3)$Sample_Well %in% p4_failed_expQC))
                                                   ))
                                                   ]

stanford_meth_rgSet_QC4
# class: RGChannelSetExtended 
# dim: 1051815 251 
# metadata(0):
#   assays(5): Green Red GreenSD RedSD NBeads
# rownames(1051815): 1600101 1600111 ... 99810990 99810992
# rowData names(0):
#   colnames(251): Plt-1-B1-DTF_het Plt-1-C1-GIST ... Plt4-G2-GIST Plt4-H2-LM
# colData names(11): Sample_Name Sample_Well ... metadata ArrayTypes
# Annotation
# array: IlluminaHumanMethylationEPIC
# annotation: ilm10b4.hg19

# stanford_meth_rgSet_QC4
# class: RGChannelSet 
# dim: 1051815 273 
# metadata(0):
#   assays(2): Green Red
# rownames(1051815): 1600101 1600111 ... 99810990 99810992
# rowData names(0):
#   colnames(273): Plt-1-B1-DTF_het Plt-1-C1-GIST ... Plt4-G2-GIST Plt4-H2-LM
# colData names(11): Sample_Name Sample_Well ... metadata ArrayTypes
# Annotation
# array: IlluminaHumanMethylationEPIC
# annotation: ilm10b4.hg19

# Save the dataset QC4 without Jurkat samples
# save(stanford_meth_rgSet_QC4, file = paste0(results_dir, "/2_stanford_meth_rgSet_BC80_logMU9_noJurkat_QC4.RData"))
# save(stanford_meth_rgSet_QC4, file = paste0(results_dir, "/2_stanford_meth_rgSetExtended_BC85_logMU9_noJurkat_QC4.RData"))
save(stanford_meth_rgSet_QC4, file = paste0(results_dir, "/2_stanford_meth_rgSetExtended_BC90_logMU9_noJurkat_QC4.RData"))

#### Calculate p values for all probes for the remaining samples
stanford_meth_pvals <- detectionP(stanford_meth_rgSet_QC4)

class(stanford_meth_pvals)
# [1] "matrix"

dim(stanford_meth_pvals)
# [1] 866091    273

# Identify genomic positions with p value > 0.01 ("failed" genomic positions
failed <- stanford_meth_pvals > 0.01

# Calculate fraction of failed positions per sample
fractions_of_failed_positions <- colMeans(failed)

# Calculate a sum of failed positions per sample
sums_of_failed_positions <- colSums(failed)

# alternative way to calculate this:
# apply(stanford_meth_pvals, 2, sum)/nrow(stanford_meth_pvals)

# Save a file with the fraction of failed positions per sample
write.table(fractions_of_failed_positions , file = paste0(qc_dir, "/pval01_fractions_of_failed_positions_per_sample_QC4.txt"),sep = "\t ", na = "NA", row.names = TRUE, col.names = TRUE)

# Save a file with the fraction of failed positions per sample
write.table(sums_of_failed_positions , file = paste0(qc_dir, "/pval01_sums_of_failed_positions_per_sample_QC4.txt"),sep = "\t ", na = "NA", row.names = TRUE, col.names = TRUE)



## Check the numbers of failed p values in samples that are technical replicates and remove the replicate that has higher number of failed p values
p1_failed_expQC <- c("B1", "C8", "F6")
p2_failed_expQC <- c("A6", "E7", "E11")
p3_failed_expQC <- c()
p4_failed_expQC <- c()

stanford_meth_rgSet_QC5 <- stanford_meth_rgSet_QC4[,
                                                   which(!(
                                                     (colData(stanford_meth_rgSet_QC4)$Sample_Plate == "Plt-1" &
                                                        (colData(stanford_meth_rgSet_QC4)$Sample_Well %in% p1_failed_expQC)) |
                                                       (colData(stanford_meth_rgSet_QC4)$Sample_Plate == "Plt2" &
                                                          (colData(stanford_meth_rgSet_QC4)$Sample_Well %in% p2_failed_expQC)) # |
                                                     #  (colData(stanford_meth_rgSet)$Sample_Plate == "Plt3" &
                                                     #  (colData(stanford_meth_rgSet)$Sample_Well %in% p3_failed_expQC)) |
                                                     # (colData(stanford_meth_rgSet)$Sample_Plate == "Plt4" &
                                                     #   (colData(stanford_meth_rgSet)$Sample_Well %in% p4_failed_expQC))
                                                   ))
                                                   ]

stanford_meth_rgSet_QC5
# class: RGChannelSetExtended 
# dim: 1051815 245 
# metadata(0):
#   assays(5): Green Red GreenSD RedSD NBeads
# rownames(1051815): 1600101 1600111 ... 99810990 99810992
# rowData names(0):
#   colnames(245): Plt-1-C1-GIST Plt-1-E1-LMS ... Plt4-G2-GIST Plt4-H2-LM
# colData names(11): Sample_Name Sample_Well ... metadata ArrayTypes
# Annotation
# array: IlluminaHumanMethylationEPIC
# annotation: ilm10b4.hg19


# stanford_meth_rgSet_QC5
# class: RGChannelSet 
# dim: 1051815 267 
# metadata(0):
#   assays(2): Green Red
# rownames(1051815): 1600101 1600111 ... 99810990 99810992
# rowData names(0):
#   colnames(267): Plt-1-C1-GIST Plt-1-D1-LMS ... Plt4-G2-GIST Plt4-H2-LM
# colData names(11): Sample_Name Sample_Well ... metadata ArrayTypes
# Annotation
# array: IlluminaHumanMethylationEPIC
# annotation: ilm10b4.hg19


# Save the dataset QC5 with technical replicate samples removed
# save(stanford_meth_rgSet_QC5, file = paste0(results_dir, "/2_stanford_meth_rgSet_BC80_logMU9_noJurkat_noReplicates_QC5.RData"))
# save(stanford_meth_rgSet_QC5, file = paste0(results_dir, "/2_stanford_meth_rgSetExtended_BC85_logMU9_noJurkat_noReplicates_QC5.RData"))
save(stanford_meth_rgSet_QC5, file = paste0(results_dir, "/2_stanford_meth_rgSetExtended_BC90_logMU9_noJurkat_noReplicates_QC5.RData"))

## Remove 17 DTF samples for the heterogeneity study 
## (multiple samples from the same patient, leave just one sample per patient with the highest number of detected probes):

p1_failed_expQC <- c("B1", "D3", "C5", "D8", "A2", "G9", "H11")
p2_failed_expQC <- c("A1", "E2", "E3", "C6", "C7", "A10", "A11")
p3_failed_expQC <- c("C9", "C10")
p4_failed_expQC <- c("B1")

stanford_meth_rgSet_QC6 <- stanford_meth_rgSet_QC5[,
                                                   which(!(
                                                     (colData(stanford_meth_rgSet_QC5)$Sample_Plate == "Plt-1" &
                                                        (colData(stanford_meth_rgSet_QC5)$Sample_Well %in% p1_failed_expQC)) |
                                                       (colData(stanford_meth_rgSet_QC5)$Sample_Plate == "Plt2" &
                                                          (colData(stanford_meth_rgSet_QC5)$Sample_Well %in% p2_failed_expQC)) |
                                                       (colData(stanford_meth_rgSet_QC5)$Sample_Plate == "Plt3" &
                                                          (colData(stanford_meth_rgSet_QC5)$Sample_Well %in% p3_failed_expQC)) |
                                                       (colData(stanford_meth_rgSet_QC5)$Sample_Plate == "Plt4" &
                                                          (colData(stanford_meth_rgSet_QC5)$Sample_Well %in% p4_failed_expQC))
                                                   ))
                                                   ]

stanford_meth_rgSet_QC6
# class: RGChannelSetExtended 
# dim: 1051815 230 
# metadata(0):
#   assays(5): Green Red GreenSD RedSD NBeads
# rownames(1051815): 1600101 1600111 ... 99810990 99810992
# rowData names(0):
#   colnames(230): Plt-1-C1-GIST Plt-1-E1-LMS ... Plt4-G2-GIST Plt4-H2-LM
# colData names(11): Sample_Name Sample_Well ... metadata ArrayTypes
# Annotation
# array: IlluminaHumanMethylationEPIC
# annotation: ilm10b4.hg19


# stanford_meth_rgSet_QC6
# class: RGChannelSet 
# dim: 1051815 251 
# metadata(0):
#   assays(2): Green Red
# rownames(1051815): 1600101 1600111 ... 99810990 99810992
# rowData names(0):
#   colnames(251): Plt-1-C1-GIST Plt-1-D1-LMS ... Plt4-G2-GIST Plt4-H2-LM
# colData names(11): Sample_Name Sample_Well ... metadata ArrayTypes
# Annotation
# array: IlluminaHumanMethylationEPIC
# annotation: ilm10b4.hg19

# Save the dataset QC6 with technical replicate samples removed
# save(stanford_meth_rgSet_QC6, file = paste0(results_dir, "/2_stanford_meth_rgSet_BC80_logMU9_noJurkat_noReplicates_DTFunique_QC6.RData"))
# save(stanford_meth_rgSet_QC6, file = paste0(results_dir, "/2_stanford_meth_rgSetExtended_BC85_logMU9_noJurkat_noReplicates_DTFunique_QC6.RData"))
save(stanford_meth_rgSet_QC6, file = paste0(results_dir, "/2_stanford_meth_rgSetExtended_BC90_logMU9_noJurkat_noReplicates_DTFunique_QC6.RData"))


## Remove 36 LM samples for the heterogeneity study 
## (multiple samples from the same patient, leave just one sample per patient with the highest number of detected probes):

p1_failed_expQC <- c("E3", "E4", "E6", "D7", "E7", "B9", "A10")
p2_failed_expQC <- c("E1", "F2", "G3", "B4", "F4", "B5", "A8", "A9", "D9", "E9", "C10", "D10", "H10", "C11", "H11", "B12", "D12", "H12")
p3_failed_expQC <- c("G1", "F4", "G4", "H4", "A5", "A7", "E7", "F7", "A8", "E8", "B12")
p4_failed_expQC <- c()

stanford_meth_rgSet_QC7 <- stanford_meth_rgSet_QC6[,
                                                   which(!(
                                                     (colData(stanford_meth_rgSet_QC6)$Sample_Plate == "Plt-1" &
                                                        (colData(stanford_meth_rgSet_QC6)$Sample_Well %in% p1_failed_expQC)) |
                                                       (colData(stanford_meth_rgSet_QC6)$Sample_Plate == "Plt2" &
                                                          (colData(stanford_meth_rgSet_QC6)$Sample_Well %in% p2_failed_expQC)) |
                                                       (colData(stanford_meth_rgSet_QC6)$Sample_Plate == "Plt3" &
                                                          (colData(stanford_meth_rgSet_QC6)$Sample_Well %in% p3_failed_expQC)) # |
                                                     # (colData(stanford_meth_rgSet_QC6)$Sample_Plate == "Plt4" &
                                                     #    (colData(stanford_meth_rgSet_QC6)$Sample_Well %in% p4_failed_expQC))
                                                   ))
                                                   ]

stanford_meth_rgSet_QC7
# class: RGChannelSetExtended 
# dim: 1051815 194 
# metadata(0):
#   assays(5): Green Red GreenSD RedSD NBeads
# rownames(1051815): 1600101 1600111 ... 99810990 99810992
# rowData names(0):
#   colnames(194): Plt-1-C1-GIST Plt-1-E1-LMS ... Plt4-G2-GIST Plt4-H2-LM
# colData names(11): Sample_Name Sample_Well ... metadata ArrayTypes
# Annotation
# array: IlluminaHumanMethylationEPIC
# annotation: ilm10b4.hg19


# stanford_meth_rgSet_QC7
# class: RGChannelSet 
# dim: 1051815 215 
# metadata(0):
#   assays(2): Green Red
# rownames(1051815): 1600101 1600111 ... 99810990 99810992
# rowData names(0):
#   colnames(215): Plt-1-C1-GIST Plt-1-D1-LMS ... Plt4-G2-GIST Plt4-H2-LM
# colData names(11): Sample_Name Sample_Well ... metadata ArrayTypes
# Annotation
# array: IlluminaHumanMethylationEPIC
# annotation: ilm10b4.hg19

# Save the dataset QC7 with multiple LM samples per patient removed
# save(stanford_meth_rgSet_QC7, file = paste0(results_dir, "/2_stanford_meth_rgSet_BC80_logMU9_noJurkat_noReplicates_DTFunique_LMunique_QC7.RData"))
# save(stanford_meth_rgSet_QC7, file = paste0(results_dir, "/2_stanford_meth_rgSetExtended_BC85_logMU9_noJurkat_noReplicates_DTFunique_LMunique_QC7.RData"))
save(stanford_meth_rgSet_QC7, file = paste0(results_dir, "/2_stanford_meth_rgSetExtended_BC90_logMU9_noJurkat_noReplicates_DTFunique_LMunique_QC7.RData"))


#####################################################################################################################################################################################
#####################################################################################################################################################################################
#####################################################################################################################################################################################
#####################################################################################################################################################################################
#####################################################################################################################################################################################


####
# Pre-Processing
####

rm(stanford_meth_rgSet)
rm(stanford_meth_rgSet_QC1)
rm(stanford_meth_rgSet_QC2)
rm(stanford_meth_rgSet_QC3)
rm(stanford_meth_rgSet_QC4)
rm(stanford_meth_rgSet_QC5)
rm(stanford_meth_rgSet_QC6)


# Remove genomic positions with probes having bead count < 3
# Seems like it needs to be done before data normalization because it is done on the RGChannelSetExtended 
# [After functional normaliztion we have a different class of the object: MethylSet]

library(wateRmelon)
# Create a matrix of bead counts with bead counts <3 represented by NA for use in the pfilter function for quality control
beadcount_matrix <- beadcount(stanford_meth_rgSet_QC7)
dim(beadcount_matrix)
# [1] 866091    194
sum(is.na(beadcount_matrix))
# [1] 279457

# Create a vector of number of samples with bead count < 3 for each probe
beadc(beadcount_matrix)

# Filter out probes with bead counts < 3 in at least one sample (>0.52% of samples, 1/194 = 0.00515)
# Set up p value filters in a way that there is no filtering done at this stage - it will be done after normalization
stanford_meth_rgSet_QC8 <- pfilter(mn = stanford_meth_rgSet_QC7, 
                                 un = stanford_meth_rgSet_QC7, 
                                 bc = beadcount_matrix, 
                                 perCount = 0.52, 
                                 perc = 10.1, 
                                 pthresh = 100, 
                                 logical.return=TRUE)
# 0 samples having 10.1 % of sites with a detection p-value greater than 0.05 were removed 
# Samples removed:  
#   50653 sites were removed as beadcount <3 in 0.52 % of samples 
# 0 sites having 100 % of samples with a detection p-value greater than 0.05 were removed 


stanford_meth_rgSet_QC8
# class: RGChannelSetExtended 
# dim: 934530 194 
# metadata(0):
#   assays(5): Green Red GreenSD RedSD NBeads
# rownames(934530): 1600101 1600111 ... 99810978 99810992
# rowData names(0):
#   colnames(194): Plt-1-C1-GIST Plt-1-E1-LMS ... Plt4-G2-GIST Plt4-H2-LM
# colData names(11): Sample_Name Sample_Well ... metadata ArrayTypes
# Annotation
# array: IlluminaHumanMethylationEPIC
# annotation: ilm10b4.hg19

save(stanford_meth_rgSet_QC8, file = paste0(results_dir, "/2_stanford_meth_rgSetExtended_BC90_logMU9_noJurkat_noReplicates_DTFunique_LMunique_nbeads3_QC8.RData"))

#### Calculate p values for all probes for the 194 samples in QC8 dataset
#### This will be used in pre-processing to remove bad quality probes
stanford_meth_pvals <- detectionP(stanford_meth_rgSet_QC8)

class(stanford_meth_pvals)
# [1] "matrix"

dim(stanford_meth_pvals)
# [1] 815438    194


## Normalization

stanford_methylSet <- preprocessFunnorm(stanford_meth_rgSet_QC8, ratioConvert = FALSE)

## Further QC
# Needed a methylSet to run this function
stanford_methylSet_QC <- minfiQC(stanford_methylSet, verbose = TRUE)

# Extract processed methylSet from returned list
stanford_methylSet <- stanford_methylSet_QC[[1]]

# Extract QC results dataframe from returned list
stanford_methylSet_QC <- stanford_methylSet_QC[[2]]; gc()

# Map to genome to create a GenomicMethylSet
stanford_gmSet <- mapToGenome(stanford_methylSet)

stanford_gmSet
# class: GenomicMethylSet 
# dim: 815231 194 
# metadata(0):
#   assays(2): Meth Unmeth
# rownames(815231): cg14817997 cg26928153 ... cg07587934 cg16855331
# rowData names(0):
#   colnames(194): Plt-1-C1-GIST Plt-1-E1-LMS ... Plt4-G2-GIST Plt4-H2-LM
# colData names(16): Sample_Name Sample_Well ... mMed uMed
# Annotation
# array: IlluminaHumanMethylationEPIC
# annotation: ilm10b4.hg19
# Preprocessing
# Method: NA
# minfi version: NA
# Manifest version: NA
# 
# save(stanford_methylSet, file = paste0(results_dir, "/3_stanford_methylSet_215samples.RData"))
# save(stanford_gmSet, file = paste0(results_dir, "/4_stanford_gmSet_215samples.RData"))
# 
# save(stanford_methylSet, file = paste0(results_dir, "/3_stanford_methylSet_Extended_209samples.RData"))
# save(stanford_gmSet, file = paste0(results_dir, "/4_stanford_gmSet_Extended_209samples.RData"))

save(stanford_methylSet, file = paste0(results_dir, "/3_stanford_methylSet_Extended_194samples.RData"))
save(stanford_gmSet, file = paste0(results_dir, "/4_stanford_gmSet_Extended_194samples.RData"))


## Visualise what the QC8 data looks like before and after normalisation #####
library(RColorBrewer)
par(mfrow=c(1,2))
densityPlot(stanford_meth_rgSet_QC8, sampGroups=colData(stanford_meth_rgSet_QC8)$Sample_Group,main="Raw", legend=FALSE, pal = brewer.pal(8, "Paired"))
legend("topleft", legend = levels(factor(colData(stanford_meth_rgSet_QC8)$Sample_Group)), 
       text.col=brewer.pal(8,"Paired"))
densityPlot(getBeta(stanford_methylSet), sampGroups=colData(stanford_meth_rgSet_QC8)$Sample_Group,
            main="Normalized", legend=FALSE, pal = brewer.pal(8, "Paired"))
legend("topleft", legend = levels(factor(colData(stanford_meth_rgSet_QC8)$Sample_Group)), 
       text.col=brewer.pal(8,"Paired"))

## Remove stanford_methylSet
rm(stanford_methylSet); gc()


## Convert to GenomicRatioSet
stanford_grSet <- ratioConvert(stanford_gmSet, what = "both")
rm(stanford_gmSet); gc()

stanford_grSet
# class: GenomicRatioSet 
# dim: 815231 194 
# metadata(0):
#   assays(3): Beta M CN
# rownames(815231): cg14817997 cg26928153 ... cg07587934 cg16855331
# rowData names(0):
#   colnames(194): Plt-1-C1-GIST Plt-1-E1-LMS ... Plt4-G2-GIST Plt4-H2-LM
# colData names(16): Sample_Name Sample_Well ... mMed uMed
# Annotation
# array: IlluminaHumanMethylationEPIC
# annotation: ilm10b4.hg19
# Preprocessing
# Method: NA
# minfi version: NA
# Manifest version: NA


## Check if sex chromosomes need to be removed
#### TODO:: determine if we want to exclude sex?
#### JP: I think yes. In the minfi 2014 paper they wrote: 
"A multidimensional scaling plot of the methylation values reveals (not surprisingly) that sex is the biggest source of variability"

#### TODO:: add condition checking if there is more than one sex in the data
#table(getSex(stanford_LMS_grSet)$predictedSex) # 32 F 45 M

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
# 18447 796784 

stanford_grSet <- stanford_grSet[keep, ]

stanford_grSet
# class: GenomicRatioSet 
# dim: 796784 194 
# metadata(0):
#   assays(3): Beta M CN
# rownames(796784): cg14817997 cg26928153 ... cg07660283 cg09226288
# rowData names(0):
#   colnames(194): Plt-1-C1-GIST Plt-1-E1-LMS ... Plt4-G2-GIST Plt4-H2-LM
# colData names(16): Sample_Name Sample_Well ... mMed uMed
# Annotation
# array: IlluminaHumanMethylationEPIC
# annotation: ilm10b4.hg19
# Preprocessing
# Method: NA
# minfi version: NA
# Manifest version: NA



## Filtering for Poor Quality Probes/Genomic positions ----

## Remove genomic positions that failed p value QC in at least 1 sample
## Calculate p values for all probes for the remaining samples


# ensure probes are in the same order in the stanford_grSet and stanford_meth_pvals objects
stanford_meth_pvals <- stanford_meth_pvals[match(featureNames(stanford_grSet),rownames(stanford_meth_pvals)),] 

# remove any probes that have failed in one or more samples
keep <- rowSums(stanford_meth_pvals < 0.01) == ncol(stanford_grSet) 
table(keep)
# keep
# FALSE   TRUE 
# 152589 644195 

stanford_grSetFlt <- stanford_grSet[keep,]

stanford_grSetFlt
# class: GenomicRatioSet 
# dim: 644195 194 
# metadata(0):
#   assays(3): Beta M CN
# rownames(644195): cg26928153 cg16269199 ... cg07660283 cg09226288
# rowData names(0):
#   colnames(194): Plt-1-C1-GIST Plt-1-E1-LMS ... Plt4-G2-GIST Plt4-H2-LM
# colData names(16): Sample_Name Sample_Well ... mMed uMed
# Annotation
# array: IlluminaHumanMethylationEPIC
# annotation: ilm10b4.hg19
# Preprocessing
# Method: NA
# minfi version: NA
# Manifest version: NA




## SNP Correction ----

# Get SNP info dataframe
stanford_SNPs <- getSnpInfo(stanford_grSetFlt)
length(grep(stanford_SNPs$CpG_rs, pattern = "rs*")) # Number of CpG SNPs
# [1] 19524
length(grep(stanford_SNPs$SBE_rs, pattern = "rs*")) # Number of SBE SNPs
# [1] 10576

# Add SNP info to grSet object
stanford_grSetFlt <- addSnpInfo(stanford_grSetFlt)

# Drop probes containing a SNP at CpG interrogation or at single NT extension
stanford_grSetFlt <- dropLociWithSnps(stanford_grSetFlt, snps=c("SBE","CpG"), maf=0)

stanford_grSetFlt
# class: GenomicRatioSet 
# dim: 624266 194 
# metadata(0):
#   assays(3): Beta M CN
# rownames(624266): cg26928153 cg16269199 ... cg07660283 cg09226288
# rowData names(6): Probe_rs Probe_maf ... SBE_rs SBE_maf
# colnames(194): Plt-1-C1-GIST Plt-1-E1-LMS ... Plt4-G2-GIST Plt4-H2-LM
# colData names(16): Sample_Name Sample_Well ... mMed uMed
# Annotation
# array: IlluminaHumanMethylationEPIC
# annotation: ilm10b4.hg19
# Preprocessing
# Method: NA
# minfi version: NA
# Manifest version: NA


## Cross-reactivity Correction

# Get list of EPIC xreactive probes
EPIC_xreactive_probes <- maxprobes::xreactive_probes()

# Filter out xreactive probes 
keep <- !(featureNames(stanford_grSetFlt) %in% EPIC_xreactive_probes)
table(keep)
# keep
# FALSE   TRUE 
# 34478 589788 

stanford_grSetFlt <- stanford_grSetFlt[keep, ]

stanford_grSetFlt
# class: GenomicRatioSet 
# dim: 589788 194 
# metadata(0):
#   assays(3): Beta M CN
# rownames(589788): cg26928153 cg16269199 ... cg06580127 cg09226288
# rowData names(6): Probe_rs Probe_maf ... SBE_rs SBE_maf
# colnames(194): Plt-1-C1-GIST Plt-1-E1-LMS ... Plt4-G2-GIST Plt4-H2-LM
# colData names(16): Sample_Name Sample_Well ... mMed uMed
# Annotation
# array: IlluminaHumanMethylationEPIC
# annotation: ilm10b4.hg19
# Preprocessing
# Method: NA
# minfi version: NA
# Manifest version: NA


save(stanford_grSetFlt, file = paste0(results_dir, "/5_stanford_grSet_filtered_preprocessed_Extended_194samples.RData"))




###############################################################################################################################################################################

## PROCESS THE 194 SAMPLE DATASET FOR EXPLORATORY ANALYSIS

########## Extract beta values for all groups together and individually
stanford_194samples_bVals <- getBeta(stanford_grSetFlt)
# Write to .csv file 
write.csv(stanford_194samples_bVals, file = paste0(results_dir, "/6_stanford_Extended_194samples_feature_beta_matrix.csv"))

## Subset the most variable beta values for unsupervised analyses (e.g. hierarchical clustering, PCA)

## Filter top 1% (5898) genomic positions with the highest variance of beta values to use for hierarchical clustering and PCA
dim(stanford_194samples_bVals)
vars=apply(stanford_194samples_bVals,1,var)
d=stanford_194samples_bVals[rev(order(vars))[1:5898],]
dim(d)

# Save matrix with top 1% most variable CpGs for the 194 sample dataset
write.csv(d, file = paste0(results_dir, "/7_stanford_Extended_194samples_feature_top1percent_beta_matrix.csv"))


######### Collapse adjacent CpG sites - Beta values
stanford_grSetFlt_Beta_collapsed <- cpgCollapse(stanford_grSetFlt, what = c('Beta'))

#Save the file with collapsed beta values:
save(stanford_grSetFlt_Beta_collapsed, file = paste0(results_dir, '/8_stanford_Extended_194samples_grSet_Beta_collapsed.RData'))

# Mapping data from collapsed sites
block_to_cpg_mapping <- stanford_grSetFlt_Beta_collapsed[[2]]$indexes

# Extract feature names
features_Beta <- rownames(stanford_194samples_bVals)

# Get the maximum number of CpG in any block
max_cpgs_per_block <- max(vapply(block_to_cpg_mapping, length, FUN.VALUE=numeric(1)))

# Rename the cpg indexes to cpg names and fill NAs up to max number of cpg
block_to_cpg_mapping_Beta <- lapply(block_to_cpg_mapping, function(block) c(features_Beta[block], rep(NA, max_cpgs_per_block - length(block))))

# Bind the list into a data.frame
block_to_cpg_mapping_Beta_df <- do.call(cbind, block_to_cpg_mapping_Beta)

# Write mapping to .csv file
write.csv(block_to_cpg_mapping_Beta_df, file= paste0(results_dir, '/9_Extended_block_to_cpg_mapping_Beta.csv'))

collapsed_Beta <- getBeta(stanford_grSetFlt_Beta_collapsed[[1]])

# Compare the dimensions of the dataset before and after collapse:
dim(stanford_grSetFlt)
# [1] 589788    194

dim(collapsed_Beta)
# [1] 322188    194

# Write to .csv file for SNF
write.csv(getBeta(stanford_grSetFlt_Beta_collapsed[[1]]), file = paste0(results_dir, "/10_Extended_stanford_194samples_v_block_Beta.csv"))



## Filter top 1% blocks (3222 blocks) with the highest variance of beta values to use for hierarchical clustering and PCA
vars=apply(collapsed_Beta,1,var)
d=collapsed_Beta[rev(order(vars))[1:3222],]
dim(d)

# Save matrix with top 1% most variable blocks
write.csv(d, file = paste0(results_dir, "/11_Extended_stanford_194samples_feature_top1percent_collapsedBeta_matrix.csv"))



###############################################################################################################################################################################

## Plot heatmaps with top 1% most variable CpGs and blocks separately

library(clustvis)

file = "7_stanford_Extended_194samples_feature_top1percent_beta_matrix.csv"
imp = importData(file)
proc = processData(imp, rowScaling = "none")
hm = generateHeatmap(proc, clustDistRows = "euclidean", clustMethodRows = "ward.D2", clustDistCols = "euclidean", clustMethodCols = "ward.D2", legendColorScheme = "Paired")
saveHeatmap(hm, file = paste0(results_dir, "/clustvisHeatmap_Extended_Euclidean_Ward2_194samples_top5percent_CpGs_01June2020.pdf"))

file = "11_Extended_stanford_194samples_feature_top1percent_collapsedBeta_matrix.csv"
imp = importData(file)
proc = processData(imp, rowScaling = "none")
hm = generateHeatmap(proc, clustDistRows = "euclidean", clustMethodRows = "ward.D2", clustDistCols = "euclidean", clustMethodCols = "ward.D2", legendColorScheme = "Paired")
saveHeatmap(hm, file = paste0(results_dir, "/clustvisHeatmap_Extended_Euclidean_Ward2_194samples_top5percent_blocks_Heatmap.pdf"))

###############################################################################################################################################################################


## Filtering to samples of interest


# All LMS samples
stanford_LMS_grSetFlt <- stanford_grSetFlt[, (colData(stanford_grSetFlt)$Sample_Group == "LMS" |
                                          colData(stanford_grSetFlt)$Sample_Group == "ULMS")]

stanford_LMS_grSetFlt
save(stanford_LMS_grSetFlt, file = paste0(results_dir, "/5a_stanford_Extended_66LMS_grSetFlt.RData" ))




############################################################################################################################################################
############################################################################################################################################################
# Differential Methylation Analysis of Regions
# - Using M value due to better distributional properties
# - Don't have a control group to compare against; can we get historic FFPE for health tissue?
############################################################################################################################################################
############################################################################################################################################################


# Extract M values (feature x sample)
stanford_LMS_M_values <- getM(stanford_LMS_grSetFlt)

# Write to .csv file for SNF
write.csv(stanford_LMS_M_values, file = paste0(results_dir, "/6a_stanford_LMS_feature_M_matrix.csv"))


# Extract beta values for all groups together and individually
stanford_LMS_bVals <- getBeta(stanford_LMS_grSetFlt)



# Write to .csv file for SNF
write.csv(stanford_LMS_bVals, file = paste0(results_dir, "/6a_stanford_LMS_feature_beta_matrix.csv"))


#########################################################################################################################################
## PROCESS SELECTED SUBSETS FOR EXPLORATORY ANALYSIS


## Subset the most variable beta values for unsupervised analyses (e.g. hierarchical clustering, PCA)

## Filter top 1000 genomic positions with the highest variance of beta values to use for hierarchical clustering and PCA
dim(stanford_LMS_bVals)
vars=apply(stanford_LMS_bVals,1,var)
d=stanford_LMS_bVals[rev(order(vars))[1:1000],]
dim(d)

# Save matrix with top 1% most variable CpGs for the 215 sample dataset
write.csv(d, file = paste0(results_dir, "/7_stanford_Extended_LMS_feature_top1000_beta_matrix.csv"))

## Filter top 1% (5898) genomic positions with the highest variance of beta values to use for hierarchical clustering and PCA
dim(stanford_LMS_bVals)
vars=apply(stanford_LMS_bVals,1,var)
d=stanford_LMS_bVals[rev(order(vars))[1:5898],]
dim(d)

# Save matrix with top 1% most variable CpGs for the 215 sample dataset
write.csv(d, file = paste0(results_dir, "/7_stanford_Extended_LMS_feature_top1percent_beta_matrix.csv"))

#########################################################################################################################################




######### Collapse adjacent CpG sites - Beta values
stanford_grSetFlt_LMS_Beta_collapsed <- cpgCollapse(stanford_LMS_grSetFlt, what = c('Beta'))

#Save the file with collapsed beta values:
save(stanford_grSetFlt_LMS_Beta_collapsed, file = paste0(results_dir, '/8_stanford_Extended_LMS_grSet_Beta_collapsed.RData'))

# Mapping data from collapsed sites
block_to_cpg_mapping <- stanford_grSetFlt_LMS_Beta_collapsed[[2]]$indexes

# Extract feature names
features_Beta <- rownames(stanford_LMS_bVals)

# Get the maximum number of CpG in any block
max_cpgs_per_block <- max(vapply(block_to_cpg_mapping, length, FUN.VALUE=numeric(1)))

# Rename the cpg indexes to cpg names and fill NAs up to max number of cpg
block_to_cpg_mapping_Beta <- lapply(block_to_cpg_mapping, function(block) c(features_Beta[block], rep(NA, max_cpgs_per_block - length(block))))

# Bind the list into a data.frame
block_to_cpg_mapping_Beta_df <- do.call(cbind, block_to_cpg_mapping_Beta)

# Write mapping to .csv file
write.csv(block_to_cpg_mapping_Beta_df, file= paste0(results_dir, '/9_Extended_block_to_cpg_mapping_LMS_Beta.csv'))

collapsed_Beta <- getBeta(stanford_grSetFlt_LMS_Beta_collapsed[[1]])

# Compare the dimensions of the dataset before and after collapse:
dim(stanford_grSetFlt)

dim(collapsed_Beta)


# Write to .csv file for SNF
write.csv(getBeta(stanford_grSetFlt_LMS_Beta_collapsed[[1]]), file = paste0(results_dir, "/10_Extended_stanford_LMS_block_Beta.csv"))





## Filter top 1% blocks (3222 blocks) with the highest variance of beta values to use for hierarchical clustering and PCA
vars=apply(collapsed_Beta,1,var)
d=collapsed_Beta[rev(order(vars))[1:3222],]
dim(d)

# Save matrix with top 1% most variable blocks
write.csv(d, file = paste0(results_dir, "/11_Extended_stanford_LMS_feature_top1percent_collapsedBeta_matrix.csv"))

vars=apply(collapsed_Beta,1,var)
d=collapsed_Beta[rev(order(vars))[1:1000],]
dim(d)

# Save matrix with top 1% most variable blocks
write.csv(d, file = paste0(results_dir, "/11_Extended_stanford_LMS_feature_top1000_collapsedBeta_matrix.csv"))

###############################################################################################################################################################################

## Plot heatmaps with top 1000 or 1% most variable CpGs and blocks separately

library(clustvis)

file = "7_stanford_Extended_LMS_feature_top1000_beta_matrix.csv"
imp = importData(file)
proc = processData(imp, rowScaling = "none")
hm = generateHeatmap(proc, clustDistRows = "euclidean", clustMethodRows = "complete", clustDistCols = "euclidean", clustMethodCols = "complete", legendColorScheme = "Paired")
saveHeatmap(hm, file = paste0(results_dir, "/clustvisHeatmap_Extended_Euclidean_complete_LMS_top1000_CpGs_01June2020.pdf"))

file = "11_Extended_stanford_LMS_feature_top1000_collapsedBeta_matrix.csv"
imp = importData(file)
proc = processData(imp, rowScaling = "none")
hm = generateHeatmap(proc, clustDistRows = "euclidean", clustMethodRows = "complete", clustDistCols = "euclidean", clustMethodCols = "complete", legendColorScheme = "Paired")
saveHeatmap(hm, file = paste0(results_dir, "/clustvisHeatmap_Extended_Euclidean_complete_LMS_top1000_blocks_Heatmap.pdf"))


file = "7_stanford_Extended_LMS_feature_top1percent_beta_matrix.csv"
imp = importData(file)
proc = processData(imp, rowScaling = "none")
hm = generateHeatmap(proc, clustDistRows = "euclidean", clustMethodRows = "complete", clustDistCols = "euclidean", clustMethodCols = "complete", legendColorScheme = "Paired")
saveHeatmap(hm, file = paste0(results_dir, "/clustvisHeatmap_Extended_Euclidean_complete_LMS_top1percent_CpGs_Heatmap.pdf"))

file = "11_Extended_stanford_LMS_feature_top1percent_collapsedBeta_matrix.csv"
imp = importData(file)
proc = processData(imp, rowScaling = "none")
hm = generateHeatmap(proc, clustDistRows = "euclidean", clustMethodRows = "complete", clustDistCols = "euclidean", clustMethodCols = "complete", legendColorScheme = "Paired")
saveHeatmap(hm, file = paste0(results_dir, "/clustvisHeatmap_Extended_Euclidean_complete_LMS_top1percent_blocks_Heatmap.pdf"))


###############################################################################################################################################################################


sessionInfo()

###
# END
###