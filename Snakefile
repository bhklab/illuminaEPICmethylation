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

import re

## TODO: Implement conda environment configuration from environment.yml
## TODO: Implement check to see if renv dependencies are available and install if not


# ---- 0. Load analysis specific configuration files and assign to local variables
configfile: 'config.yaml'

# Unpack config dictionary to local environment
plate_data_dirs = config['plate_data_dirs']
plate_labels = config['plate_labels']

nthread = config['nthread']

analysis_name = config['analysis_name']

detection_pvalue = config['detection_pvalue']

bisulphite_conversion_rate = config['bisulphite_conversion_rate']


# ---- 1. Build minfi RGSet from raw plate data and labels

rule build_rgset_from_plate_data:
    input:
        plates=expand('rawdata/{plate_dirs}', plate_dirs=plate_data_dirs),
        labels=expand('metadata/{labels}', labels=plate_labels)
    output:
        rgset=f'procdata/1.{analysis_name}.RGChannelSet.qs'
    threads: nthread
    script:
        "scripts/1_buildRGsetFromPlateData.R"


# ---- 2. Filter out samples with less than 90% CpG probes detected

qc_path1 = f'qc/2.{analysis_name}.RGChannelSet'
qc_steps1 = [f'detectionP{detection_pvalue}']
qc_step1 = qc_steps1[0]

rule drop_samples_less_than_90pct_probes_detected:
    input:
        rgset=f'procdata/1.{analysis_name}.RGChannelSet.qs'
    params:
        detection_pvalue=detection_pvalue
    output:
        detectionPvalues=f'{qc_path1}.detection_pvalues.csv',
        probeQC=f'{qc_path1}.num_probes_with_proportion_failed_samples_p{detection_pvalue}.csv',
        sampleQC=f'{qc_path1}.probes_failed_per_sample_p{detection_pvalue}.csv',
        rgset_filtered=f'procdata/2.{analysis_name}.RGChannelSet.{qc_step1}.qs'
    script:
        "scripts/2_dropSamplesLessThan90PctProbesDetected.R"


# ---- 3. Filter out samples with bisulphite conversion rate lower than specified in 'config.yml'

qc_path2 = f'qc/3.{analysis_name}.RGChannelSet'
qc_step2_list = [*qc_steps1, f'bisulphite{bisulphite_conversion_rate}']
qc_step2 = '.'.join(qc_step2_list)

rule drop_samples_with_bisulphite_conversion_less_than_cutoff:
    input:
        rgset=f'procdata/2.{analysis_name}.RGChannelSet.{qc_step1}.qs'
    params:
        bisulphite_conversion_rate=bisulphite_conversion_rate
    output:
        bisulphite_qc=f'{qc_path2}.bisulphite_conversions.csv',
        rgset_filtered=f'procdata/3.{analysis_name}.RGChannelSet.{qc_step2}.qs'
    script:
        "scripts/3_dropSamplesBisulfiteConversionLessThanCutoff.R"


# ---- 4. Filter out samples failing one or more manual QC steps specified in 'config.yml'

# Read in relevant configuration
manual_qc_steps = config['manual_qc_steps']

# Build manual qc output file names
manual_qc_step_names = list(manual_qc_steps.keys())

manual_qc_file_name = f'procdata/4.{analysis_name}.RGChannelSet'
manual_qc_output_file_names = []
qc_step3_list = qc_step2_list

for step in manual_qc_step_names:
    qc_step3_list = [*qc_step3_list, step]
    manual_qc_current_step = '.'.join(qc_step3_list)
    manual_qc_output_file_names = [*manual_qc_output_file_names, 
                                  f'{manual_qc_file_name}.{manual_qc_current_step}.qs']

qc_step_num = 1 + len(manual_qc_step_names) 

rule drop_samples_failed_manual_qc:
    input:
        rgset=f'procdata/3.{analysis_name}.RGChannelSet.{qc_step2}.qs'
    params:
        manual_qc_steps=manual_qc_steps
    output:
        file_names=manual_qc_output_file_names
    script:
        "scripts/4_dropSamplesFailedManualQC.R"


# ---- 5. Filter out probes with less than 3 beads

qc_steps_final = [*qc_step3_list, 'probesGt3Beads']
final_qc_step = '.'.join(qc_steps_final)

rule drop_probes_with_less_than_three_beads:
    input:
        rgset=f'procdata/4.{analysis_name}.RGChannelSet.{manual_qc_current_step}.qs'
    output:
        bead_counts=f'qc/5.{analysis_name}.RGChannelSet.{manual_qc_current_step}.bead_counts.csv',
        rgset_filtered=f'procdata/5.{analysis_name}.RGChannelSet.{final_qc_step}.qs'
    shell:
        """
        Rscript scripts/5_dropProbesWithLessThan3Beads.R \
            -r {input.rgset} \
            -b {output.bead_counts} \
            -o {output.rgset_filtered}
        """


# ---- 6. Preprocess RGChannelSet to MethylSet

# Read in relevant configuration
preprocess_methods = config['preprocess_methods']

rule preprocess_to_methylset_and_qc:
    input:
        rgset=f'procdata/5.{analysis_name}.RGChannelSet.{final_qc_step}.qs'
    params:
        preprocess_methods=preprocess_methods
    output:
        methylsets=expand('procdata/6.{analysis_name}.MethylSet.{preprocess_method}.qs',
            analysis_name=analysis_name, preprocess_method=preprocess_methods),
        qc_reports=expand('qc/6.{analysis_name}.MethylSet.{preprocess_method}.qc_report.csv',
            analysis_name=analysis_name, preprocess_method=preprocess_methods)
    threads: nthread
    script:
        'scripts/6_preprocessToMethylSetAndQC.R'



# ---- 7. Visual the Raw RGSet vs Each Preprocessing Method with QC Metrics

preprocess_string = '_vs_'.join(preprocess_methods_split)
comparisons = f'rgSet_vs_{preprocess_string}'

rule density_plot_preprocessed_vs_rgset:
    input:
        methylsets=expand('procdata/6.{analysis_name}.MethylSet.{preprocess_method}.qs',
            analysis_name=analysis_name, preprocess_method=preprocess_methods_split),
        rgset=f'procdata/5.{analysis_name}.RGChannelSet.{final_qc_step}.qs'
    output:
        plots=f'qc/7.{analysis_name}.{comparisons}.density_plots.pdf',
        qc_report=f'qc/7.{analysis_name}.{comparisons}.stats_table.csv'
    threads: nthread
    shell:
        """
        Rscript scripts/7_densityPlotPreprocessedVsRGSet.R \
            -m '{input.methylsets}' \
            -r {input.rgset} \
            -o '{output.plots}' \
            -q '{output.qc_report}' \
            -t '{preprocess_methods}'
        """


# ---- 8. Filter out one or more manual QC steps after prepreprocessing

# Read in configuration for this step
selected_preprocess_method = config['selected_preprocess_method']
manual_qc2_steps = list(config['manual_qc2_steps'].keys())

# Build outpuit file names for iteratively strict qc steps
manual_qc2_file_name = f'procdata/8.{analysis_name}.MethylSet.{selected_preprocess_method}'
manual_qc2_output_file_names = []
qc2_steps = []

for step in manual_qc2_steps:
    qc2_steps = [*qc2_steps, step]
    manual_qc2_current_step = '.'.join(qc2_steps)
    manual_qc2_output_file_names = [*manual_qc2_output_file_names, 
                                    f'{manual_qc2_file_name}.{manual_qc2_current_step}.qs']


rule drop_samples_failed_manual_qc2:
    input:
        methylset=f'procdata/6.{analysis_name}.MethylSet.{selected_preprocess_method}.qs',
        qc_report=f'qc/7.{analysis_name}.{comparisons}.stats_table.csv'
    params:
        preproc_method=selected_preprocess_method,
        qc_criteria=config['manual_qc2_steps']
    output:
        manual_qc2_output_file_names
    threads: nthread
    script:
        "scripts/8_dropSamplesFailedManualQC2.R"


# ---- 9. Convert GenomicMethylSet to GenomicRatioSet and add annotations

## TODO: Implement support for multiple QC2 step selection

selected_qc2_step = config['selected_qc2_step']

rule convert_gmset_to_grset_and_drop_sex_chromosomes:
    input:
        methylset=f'procdata/8.{analysis_name}.MethylSet.{selected_preprocess_method}.{selected_qc2_step}.qs'
    output:
        genomicmethylset=f'procdata/9.{analysis_name}.{selected_preprocess_method}.GenomicMethylSet.qs',
        ratioset=f'procdata/9.{analysis_name}.{selected_preprocess_method}.GenomicRatioSet.drop_sex_chr.qs'
    threads: nthread
    shell:
        """
        Rscript scripts/9_convertGMSetToGRSetAndDropSexChromosomes.R \
            -i {input.methylset} \
            -a {output.genomicmethylset} \
            -o {output.ratioset}
        """


# ---- 10. Filter poor quality probes
rule filter_grset_poor_quality_probes:
    input:
        grset=f'procdata/9.{analysis_name}.{selected_preprocess_method}.GenomicRatioSet.drop_sex_chr.qs',
        pvalues=f'qc/2.{analysis_name}.RGChannelSet.detection_pvalues.csv'
    output:
        filtered_grset=f'procdata/10.{analysis_name}.{selected_preprocess_method}.GenomicRatioSet.drop_sex_chr.filter_probes.qs'
    threads: nthread
    shell:
        """
        Rscript scripts/10_filterGRSetPoorQualityProbes.R \
            -g {input.grset} \
            -p {input.pvalues} \
            -o {output.filtered_grset}
        """


# ---- 11. SNP correction

rule correct_grset_for_snps:
    input:
        grset=f'procdata/10.{analysis_name}.{selected_preprocess_method}.GenomicRatioSet.drop_sex_chr.filter_probes.qs'
    output:
        drop_snps_grset=f'procdata/11.{analysis_name}.{selected_preprocess_method}.GenomicRatioSet.drop_sex_chr.filter_probes.drop_snps.qs'
    threads: nthread
    shell:
        """
        Rscript scripts/11_correctGRSetForSNPs.R \
            -g {input.grset} \
            -o {output.drop_snps_grset}
        """


# ---- 12. Correct for cross-reactive probes

rule correct_grset_for_crossreactive_probes:
    input:
        grset=f'procdata/11.{analysis_name}.{selected_preprocess_method}.GenomicRatioSet.drop_sex_chr.filter_probes.drop_snps.qs'
    output:
        drop_xreactive_grset=f'procdata/12.{analysis_name}.{selected_preprocess_method}.GenomicRatioSet.drop_sex_chr.filter_probes.drop_snps.correct_xreactive.qs'
    threads: nthread
    shell:
        """
        Rscript scripts/12_correctGRSetForCrossReactiveProbes.R \
            -g {input.grset} \
            -o {output.drop_xreactive_grset}
        """

## TODO:: Iteratively build the file labels so there aren't so many parameters being passed to input/output

# ---- 13. Subset grSet by cancer type

cancer_types = config['cancer_types']

rule subset_grset_by_cancer_types:
    input:
        grset=f'procdata/12.{analysis_name}.{selected_preprocess_method}.GenomicRatioSet.drop_sex_chr.filter_probes.drop_snps.correct_xreactive.qs'
    output:
        grset=f'results/13.{analysis_name}.{selected_preprocess_method}.all_types.GenomicRatioSet.qs',
        grsets=expand('results/13.{analysis_name}.{selected_preprocess_method}.{cancer_type}.GenomicRatioSet.qs', 
                      analysis_name=analysis_name, cancer_type=cancer_types, selected_preprocess_method=selected_preprocess_method)
    threads: nthread
    shell:
        """
        Rscript scripts/13_subsetSamplesByCancerType.R \
            -g {input.grset} \
            -s '{cancer_types}' \
            -o '{output.grset} {output.grsets}'
        """


cancer_types_all = ['all_types', *cancer_types] 

# ---- 14. Extract M and Beta values from each GenomicRatioSet for CpGs

rule extract_m_and_beta_values:
    input:
        grsets=expand('results/13.{analysis_name}.{selected_preprocess_method}.{cancer_type}.GenomicRatioSet.qs', 
                      analysis_name=analysis_name, cancer_type=cancer_types_all, selected_preprocess_method=selected_preprocess_method)
    output:
        m_values=expand('results/14.{analysis_name}.{selected_preprocess_method}.{cancer_type}.m_values.csv',
                         analysis_name=analysis_name, cancer_type=cancer_types_all, selected_preprocess_method=selected_preprocess_method),
        beta_values=expand('results/14.{analysis_name}.{selected_preprocess_method}.{cancer_type}.beta_values.csv',
                           analysis_name=analysis_name, cancer_type=cancer_types_all, selected_preprocess_method=selected_preprocess_method)
    threads: nthread
    shell:
        """
        Rscript scripts/14_extractCpGMandBetaValues.R \
            -g '{input.grsets}' \
            -b '{output.beta_values}' \
            -m '{output.m_values}'
        """


# --- 15. Collapse Adjacent CpG sites into methyalted regions for M and Beta values
rule collapse_grset_cpgs_to_methylated_regions:
    input:
        grsets=expand('results/13.{analysis_name}.{selected_preprocess_method}.{cancer_type}.GenomicRatioSet.qs', 
                        analysis_name=analysis_name, cancer_type=cancer_types_all, selected_preprocess_method=selected_preprocess_method)
    output:
        m_region_grsets=expand('results/15.{analysis_name}.{selected_preprocess_method}.{cancer_type}.GenomicRatioSet.m_values.regions.qs',
                                 analysis_name=analysis_name, cancer_type=cancer_types_all, selected_preprocess_method=selected_preprocess_method),
        beta_region_grsets=expand('results/15.{analysis_name}.{selected_preprocess_method}.{cancer_type}.GenomicRatioSet.beta_values.regions.qs',
                                  analysis_name=analysis_name, cancer_type=cancer_types_all, selected_preprocess_method=selected_preprocess_method)
    threads: nthread
    shell:
        """
        Rscript scripts/15_collapseGRSetCpGsToMethylatedRegions.R \
            -g '{input.grsets}' \
            -b '{output.beta_region_grsets}' \
            -m '{output.m_region_grsets}'
        """

# ---- 16. Build methylated region to CpG mappings and write region M and Beta values to disk

methylation_values = ['m_values', 'beta_values']

rule build_region_mappings_for_beta_and_m_values:
    input:
        grsets=expand('results/15.{analysis_name}.{selected_preprocess_method}.{cancer_type}.GenomicRatioSet.{methylation_values}.regions.qs',
            analysis_name=analysis_name, cancer_type=cancer_types_all, methylation_values=methylation_values, selected_preprocess_method=selected_preprocess_method),
        methylation_data=expand('results/14.{analysis_name}.{selected_preprocess_method}.{cancer_type}.{methylation_values}.csv',
            analysis_name=analysis_name, cancer_type=cancer_types_all, methylation_values=methylation_values, selected_preprocess_method=selected_preprocess_method)
    output:
        mappings=expand('results/16.{analysis_name}.{selected_preprocess_method}.{cancer_type}.{methylation_value}.region_to_cpg_mappings.csv', 
                        analysis_name=analysis_name, cancer_type=cancer_types_all, methylation_value=methylation_values, selected_preprocess_method=selected_preprocess_method),
        methyl_values=expand('results/16.{analysis_name}.{selected_preprocess_method}.{cancer_type}.{methylation_value}.regions.csv', 
                             analysis_name=analysis_name, cancer_type=cancer_types_all, methylation_value=methylation_values, selected_preprocess_method=selected_preprocess_method)
    threads: nthread
    shell:
        """
        Rscript scripts/16_buildRegionMappingsForBetaAndMValues.R \
            -g '{input.grsets}' \
            -d '{input.methylation_data}' \
            -M '{output.mappings}' \
            -V '{output.methyl_values}'
        """