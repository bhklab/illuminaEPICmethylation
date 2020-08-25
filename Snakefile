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

failed_qc1=config['failed_qc1']
failed_qc2=config['failed_qc2']
failed_qc3=config['failed_qc3']


# ---- 1. Build minfi RGSet from raw plate data and labels

rule build_rgset_from_plate_data:
    input:
        plates=expand('rawdata/{plate_dirs}', plate_dirs=plate_data_dirs),
        labels=expand('metadata/{labels}', labels=plate_labels)
    output:
        expand('procdata/1_{analysis_name}_RGSet_raw.qs', analysis_name=analysis_name)
    shell:
        """
        Rscript scripts/buildRGsetFromPlateData.R \
            -p '{input.plates}' \
            -l '{input.labels}' \
            -o '{output}' \
            -n {nthread}
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
        sample_qc=expand('qc/2_{analysis_name}_num_samples_with_proportion_failed_probes.csv', analysis_name=analysis_name),
        qc_report=expand('qc/2_{analysis_name}_minfi_qc_report.pdf', analysis_name=analysis_name)
    shell:
        """
        Rscript scripts/generateMicroarrayQCReport.R \
            -i {input.rgset} -d {output.detection_pvals} \
            -p {output.probe_qc} \
            -s {output.sample_qc} \
            -r {output.qc_report}
        """


# ---- 3. Preprocess microarray intensities to methylation values

rule convert_rgset_to_methylset_for_qc:
    input:
        rgset=expand('procdata/1_{analysis_name}_RGSet_raw.qs', analysis_name=analysis_name)
    output:
        qc_figures=expand('qc/methylSet/3_{analysis_name}_methylset_preproc_qc_plots.pdf', analysis_name=analysis_name),
        methylset=expand('procdata/3_{analysis_name}_methylset.qs', analysis_name=analysis_name)
    shell:
        """
        Rscript scripts/preprocessRGSetToMethylSet.R \
            -i {input.rgset} \
            -f {output.qc_figures} \
            -o {output.methylset}
        """


## FIXME: How can I reuse a rule with different outputs?
## TODO: Determine if we can automate any manual QC steps

# ---- 4. Pre-normalization manual QC1

## Outlier samples with low signal from methylated and unmethylated channels, 
## low values of control probes and poor density plots: Removed 15 samples that 
## failed at least two different QC metrics

rule remove_failed_qc1_and_generate_report:
    input:
        rgset=expand('procdata/1_{analysis_name}_RGSet_raw.qs', analysis_name=analysis_name)
    output:
        rgset_qc=expand('procdata/4_{analysis_name}_RGSet_qc1.qs', analysis_name=analysis_name),
        detection_pvals=expand('qc/prenormalization/4_{analysis_name}_qc1_detection_pvals.csv', analysis_name=analysis_name),
        probe_qc=expand('qc/prenormalization/4_{analysis_name}__qc1_probes_failed_per_sample_p0.01.csv', analysis_name=analysis_name),
        sample_qc=expand('qc/prenormalization/4_{analysis_name}_qc1_num_samples_with_proportion_failed_probes.csv', analysis_name=analysis_name)
    shell:
        """
        Rscript scripts/removeFailedQCandGenerateReport.R \
            -i {input.rgset} \
            -f '{failed_qc1}' \
            -d {output.detection_pvals} \
            -p {output.probe_qc} \
            -s {output.sample_qc} \
            -o {output.rgset_qc}
        """

# ---- 5. Pre-normalization manual QC2

rule remove_failed_qc2_and_generate_report:
    input:
        rgset=expand('procdata/4_{analysis_name}_RGSet_qc1.qs', analysis_name=analysis_name)
    output:
        rgset_qc2=expand('procdata/5_{analysis_name}_RGSet_qc2.qs', analysis_name=analysis_name),
        detection_pvals=expand('qc/prenormalization/5_{analysis_name}_qc2_detection_pvals.csv', analysis_name=analysis_name),
        probe_qc=expand('qc/prenormalization/5_{analysis_name}__qc2_probes_failed_per_sample_p0.01.csv', analysis_name=analysis_name),
        sample_qc=expand('qc/prenormalization/5_{analysis_name}_qc2_num_samples_with_proportion_failed_probes.csv', analysis_name=analysis_name)    
    shell:
        """
        Rscript scripts/removeFailedQCandGenerateReport.R \
            -i {input.rgset} \
            -f '{failed_qc2}' \
            -d {output.detection_pvals} \
            -p {output.probe_qc} \
            -s {output.sample_qc} \
            -o {output.rgset_qc2}
        """

# ---- 6. Pre-normalization manual QC3

rule remove_failed_qc3_and_generate_report:
    input:
        rgset=expand('procdata/5_{analysis_name}_RGSet_qc2.qs', analysis_name=analysis_name)
    output:
        rgset_qc3=expand('procdata/6_{analysis_name}_RGSet_qc3.qs', analysis_name=analysis_name),
        detection_pvals=expand('qc/prenormalization/6_{analysis_name}_detection_pvals.csv', analysis_name=analysis_name),
        probe_qc=expand('qc/prenormalization/6_{analysis_name}__qc3_probes_failed_per_sample_p0.01.csv', analysis_name=analysis_name),
        sample_qc=expand('qc/prenormalization/6_{analysis_name}_qc3_num_samples_with_proportion_failed_probes.csv', analysis_name=analysis_name)
    shell:
        """
        Rscript scripts/removeFailedQCandGenerateReport.R \
            -i {input.rgset} \
            -f '{failed_qc3}' \
            -d {output.detection_pvals} \
            -p {output.probe_qc} \
            -s {output.sample_qc} \
            -o {output.rgset_qc3}
        """

# ---- 7. Functional normalize RGSet

rule functional_normalize_rgset_to_methylset:
    input:
        rgset=f'procdata/6_{analysis_name}_RGSet_qc3.qs'
    output:
        methylset=f'procdata/7_{analysis_name}_methylset_raw.qs',
        qc_report=f'qc/normalized/7_{analysis_name}_methylset_qc1.csv',
        plot=f'qc/normalized/7_{analysis_name}_normalized vs unnormalized_distributions.pdf'
    shell:
        """
        Rscript scripts/functionalNormalizeAndQC.R \
            -i {input.rgset} \
            -o {output.methylset} \
            -r {output.qc_report}
        """

# ---- 8. Visualize normalied data vs QC2 and QC3 unnormalized data

rule plot_normalized_vs_qc2_and_qc3:
    input:
        rgset_qc2=f'procdata/5_{analysis_name}_RGSet_qc2.qs',
        rgset_qc3=f'procdata/6_{analysis_name}_RGSet_qc3.qs',
        normalized=f'procdata/7_{analysis_name}_methylset_raw.qs'
    output:
        plot1=f'qc/normalized/8_{analysis_name}_normalied_vs_unnormalized_QC2.pdf',
        plot2=f'qc/normalized/8_{analysis_name}_normalied_vs_unnormalized_QC3.pdf'
    shell:
        """
        Rscript scripts/plotNormalizedVsQc2AndQc3.R \
            -q '{input.rgset_qc2} {input.rgset_qc3}' \
            -n {input.normalized} \
            -o '{output.plot1} {output.plot2}'
        """

# ---- 9. Convert GenomicMethylSet to GenomicRatioSet and add annotations

rule convert_gmset_to_grset_and_annotate:
    input:
    output:
    shell:
        """
        Rscript convertGMSetToGRSetAndAnnotate.R \
            - i \
            - o 
        """


# ---- 10. Filter poor quality probes
rule filter_grset_poor_quality_probes:
    input:
    output:
    shell:
        """
        Rscript filterGRSetPoorQualityProbes \
            -i 
        """


# ---- 11. SNP correction

rule correct_grset_for_snps:
    input:
    output:
    shell:
        """
        Rscript correctGRSetForSNPs.R \
            -i  
        """

# ---- 12. Correct for cross-reactive probes

rule correct_grset_for_crossreactive_probes:
    input:
    output:
    shell:
        """
        Rscript correctGRSetForCrossReactiveProbesS.R \
            -i 
        """
