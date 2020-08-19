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
        Rscript scripts/generateMicroarrayQCReport.R -i {input.rgset} \\
            -d {output.detection_pvals} -p {output.probe_qc} \\ 
            -s {output.sample_qc} -r {output.qc_report}
        """

# ---- 3. 