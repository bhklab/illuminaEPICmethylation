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

import re

## TODO: Implement conda environment configuration from environment.yml
## TODO: Implement check to see if packrat dependencies are available and install if not


# ---- 0. Load analysis specific configuration files and assign to local variables
configfile: 'config.yaml'

plate_data_dirs = config['plate_data_dirs']
plate_labels = config['plate_labels']
nthread = config['nthread']
analysis_name = config['analysis_name']


# ---- 1. Build minfi RGSet from raw plate data and labels

rule build_rgset_from_plate_data:
    input:
        plates=expand('rawdata/{plate_dirs}', plate_dirs=plate_data_dirs),
        labels=expand('metadata/{labels}', labels=plate_labels)
    output:
        expand('procdata/1_{analysis_name}_RGSet_raw.qs', analysis_name=analysis_name)
    shell:
        """
        Rscript scripts/buildRGsetFromPlateData.R -p '{input.plates}' -l '{input.labels}' -o '{output}' -n {nthread}
        """


# ---- 2. Generate microarray QC report

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

rule generate_microarray_qc_report:
    input:
        rgset=expand('procdata/1_{analysis_name}_RGSet_raw.qs', analysis_name=analysis_name)
    output:
        detection_pvals=expand('qc/2_{analysis_name}_detection_pvals.csv', analysis_name=analysis_name),
        probe_qc=expand('qc/2_{analysis_name}_probes_failed_per_sample_p0.01.csv', analysis_name=analysis_name),
        sample_qc=expand('qc/2_{analysis_name}_num_samples_failed_at_pval_cutoff.csv', analysis_name=analysis_name),
        qc_report=expand('qc/2_{analysis_name}_minfi_qc_report.pdf', analysis_name=analysis_name)
    shell:
        """
        Rscript scripts/generateMicroarrayQCReport.R -i {input.rgset} -d {output.detection_pvals} -p {output.probe_qc} -s {output.sample_qc} -r {output.qc_report}
        """


# ---- 3. Preprocess microarray intensities to methylation values

rule preprocess_rgset_to_methylset:
    input:
        rgset=expand('procdata/1_{analysis_name}_RGSet_raw.qs', analysis_name=analysis_name)
    output:
        qc_figures=expand('qc/methylSet/3_{analysis_name}_methylset_preproc_qc_plots.pdf', analysis_name=analysis_name),
        methylset=expand('procdata/3_{analysis_name}_methylset.qs', analysis_name=analysis_name)
    shell:
        """
        Rscript scripts/preprocessRGSetToMethylSet.R -i {input.rgset} -f {output.qc_figures} -o {output.methylset}
        """