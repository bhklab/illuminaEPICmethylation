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

manual_qc_steps = config['manual_qc_steps']

preprocess_methods = config['preprocess_methods']

manual_qc2_steps = config['manual_qc2_steps']

cancer_types=config['cancer_types']


# ---- 1. Build minfi RGSet from raw plate data and labels

rule build_rgset_from_plate_data:
    input:
        plates=expand('rawdata/{plate_dirs}', plate_dirs=plate_data_dirs),
        labels=expand('metadata/{labels}', labels=plate_labels)
    output:
        f'procdata/1.{analysis_name}.RGChannelSet.qs'
    threads: nthread
    shell:
        """
        Rscript scripts/1_buildRGsetFromPlateData.R \
            -p '{input.plates}' \
            -l '{input.labels}' \
            -o '{output}' \
            -n {nthread}
        """


# ---- 2. Filter out samples with less than 90% CpG probes detected

qc_path1 = f'qc/2.{analysis_name}.RGChannelSet'
qc_steps1 = [f'detectionP{detection_pvalue}']
qc_step1 = qc_steps1[0]

rule drop_samples_less_than_90pct_probes_detected:
    input:
        rgset=f'procdata/1.{analysis_name}.RGChannelSet.qs'
    output:
        detectionPvalues=f'{qc_path1}.detection_pvalues.csv',
        probeQC=f'{qc_path1}.num_probes_with_proportion_failed_samples_p{detection_pvalue}.csv',
        sampleQC=f'{qc_path1}.probes_failed_per_sample_p{detection_pvalue}.csv',
        rgset_filtered=f'procdata/2.{analysis_name}.RGChannelSet.{qc_step1}.qs'
    shell:
        """
        Rscript scripts/2_dropSamplesLessThan90PctProbesDetected.R \
            -r {input.rgset} \
            -p {detection_pvalue} \
            -P {qc_path1} \
            -o {output.rgset_filtered}
        """


# ---- 3. Filter out samples with bisulphite conversion rate lower than specified in 'config.yml'

qc_path2 = f'qc/3.{analysis_name}.RGChannelSet'
qc_steps2 = [*qc_steps1, f'bisulphite{bisulphite_conversion_rate}']
qc_step2 = '.'.join(qc_steps2)

rule drop_samples_with_bisulphite_conversion_less_than_cutoff:
    input:
        rgset=f'procdata/2.{analysis_name}.RGChannelSet.{qc_step1}.qs'
    output:
        bisulphite_conversion=f'{qc_path2}.bisulphite_conversions.csv',
        rgset_filtered=f'procdata/3.{analysis_name}.RGChannelSet.{qc_step2}.qs'
    shell:
        """
        Rscript scripts/3_dropSamplesBisulfiteConversionLessThanCutoff.R\
            -r {input.rgset} \
            -c {bisulphite_conversion_rate} \
            -P {qc_path2} \
            -o {output.rgset_filtered}
        """


# ---- 4. Filter out samples failing one or more manual QC steps specified in 'config.yml'

# Build manual qc output file names
manual_qc_step_names = [re.sub(':.*$', '', step) for step in manual_qc_steps]

manual_qc_file_name = f'procdata/4.{analysis_name}.RGChannelSet'
manual_qc_file_names = []
qc_steps3 = qc_steps2

for step in manual_qc_step_names:
    qc_steps3 = [*qc_steps3, step]
    manual_qc_step = '.'.join(qc_steps3)
    manual_qc_file_names = [*manual_qc_file_names, f'{manual_qc_file_name}.{manual_qc_step}.qs']

qc_step_num = 1 + len(manual_qc_step_names) 

rule drop_samples_failed_manual_qc:
    input:
        rgset=f'procdata/3.{analysis_name}.RGChannelSet.{qc_step2}.qs'
    output:
        manual_qc_file_names
    shell:
        """
        Rscript scripts/4_dropSamplesFailedManualQC.R \
            -r {input.rgset} \
            -q '{manual_qc_steps}' \
            -o '{output}'
        """


# ---- 5. Filter out probes with less than 3 beads

qc_steps_final = [*qc_steps3, 'probesGt3Beads']
final_qc_step = '.'.join(qc_steps_final)

rule drop_probes_with_less_than_three_beads:
    input:
        rgset=f'procdata/4.{analysis_name}.RGChannelSet.{manual_qc_step}.qs'
    output:
        bead_counts=f'qc/5.{analysis_name}.RGChannelSet.{manual_qc_step}.bead_counts.csv',
        rgset_filtered=f'procdata/5.{analysis_name}.RGChannelSet.{final_qc_step}.qs'
    shell:
        """
        Rscript scripts/5_dropProbesWithLessThan3Beads.R \
            -r {input.rgset} \
            -b {output.bead_counts} \
            -o {output.rgset_filtered}
        """


# ---- 6. Functional normalize RGChannelSet to MethylSet
preprocess_methods = preprocess_methods.split(',')

rule preprocess_to_methylset_and_qc:
    input:
        rgset=f'procdata/5.{analysis_name}.RGChannelSet.{final_qc_step}.qs'
    output:
        methylsets=expand('procdata/6.{analysis_name}.MethylSet.{preprocess_method}.qs',
            analysis_name=analysis_name, preprocess_method=preprocess_methods),
        qc_reports=expand('qc/6.{analysis_name}.MethylSet.{preprocess_method}.qc_report.csv',
            analysis_name=analysis_name, preprocess_method=preprocess_methods)
    threads: nthread
    shell:
        """
        Rscript scripts/6_preprocessToMethylSetAndQC.R \
            -i {input.rgset} \
            -m '{preprocess_methods}' \
            -o '{output.methylsets}' \
            -r '{output.qc_reports}'
        """


# ---- 7. Visualize normalized vs unnormalized beta value distribution

rule plot_normalized_vs_two_previous_qc_steps:
    input:
        rgset=f'procdata/4.{analysis_name}.RGChannelSet.{final_qc_step}.qs',
        normalized=expand('procdata/5.{analysis_name}.MethylSet.{preprocess_method}.qs',
            analysis_name=analysis_name, preprocess_method=preprocess_methods)
    output:
        plot=f'qc/6.{analysis_name}.RGChannelSet.{final_qc_step}.vs_normalized_plot.pdf',
    threads: nthread
    shell:
        """
        Rscript scripts/6_normalizedVs2PreviousQCSteps.R \
            -q {input.rgset} \
            -n {input.normalized} \
            -o {output.plot}
        """


# ---- 9. Convert GenomicMethylSet to GenomicRatioSet and add annotations

rule convert_gmset_to_grset_and_drop_sex_chromosomes:
    input:
        methylset=f'procdata/7.{analysis_name}.MethylSet.funnorm.qs'
    output:
        genomicmethylset=f'procdata/9.{analysis_name}.GenomicMethylSet.qs',
        ratioset=f'procdata/9.{analysis_name}.GenomicRatioSet.drop_sex_chr.qs'
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
        grset=f'procdata/9.{analysis_name}.GenomicRatioSet.drop_sex_chr.qs',
        pvalues=f'qc/6.{analysis_name}.RGChannelSet.qc3.detection_pvals.csv'
    output:
        filtered_grset=f'procdata/10.{analysis_name}.GenomicRatioSet.drop_sex_chr.filter_probes.qs'
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
        grset=f'procdata/10.{analysis_name}.GenomicRatioSet.drop_sex_chr.filter_probes.qs'
    output:
        drop_snps_grset=f'procdata/11.{analysis_name}.GenomicRatioSet.drop_sex_chr.filter_probes.drop_snps.qs'
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
        grset=f'procdata/11.{analysis_name}.GenomicRatioSet.drop_sex_chr.filter_probes.drop_snps.qs'
    output:
        drop_xreactive_grset=f'procdata/12.{analysis_name}.GenomicRatioSet.drop_sex_chr.filter_probes.drop_snps.correct_xreactive.qs'
    threads: nthread
    shell:
        """
        Rscript scripts/12_correctGRSetForCrossReactiveProbes.R \
            -g {input.grset} \
            -o {output.drop_xreactive_grset}
        """


# ---- 13. Subset grSet by cancer type

rule subset_grset_by_cancer_types:
    input:
        grset=f'procdata/12.{analysis_name}.GenomicRatioSet.drop_sex_chr.filter_probes.drop_snps.correct_xreactive.qs'
    output:
        grset=f'results/13.{analysis_name}.all_types.GenomicRatioSet.qs',
        grsets=expand('results/13.{analysis_name}.{cancer_type}.GenomicRatioSet.qs', 
                      analysis_name=analysis_name, cancer_type=cancer_types)
    threads: nthread
    shell:
        """
        Rscript scripts/13_subsetSamplesByCancerType.R \
            -g {input.grset} \
            -s '{cancer_types}' \
            -o '{output.grset} {output.grsets}'
        """


cancer_types = ['all_types', *cancer_types] 

# ---- 14. Extract M and Beta values from each GenomicRatioSet for CpGs

rule extract_m_and_beta_values:
    input:
        grsets=expand('results/13.{analysis_name}.{cancer_type}.GenomicRatioSet.qs', 
                      analysis_name=analysis_name, cancer_type=cancer_types)
    output:
        m_values=expand('results/14.{analysis_name}.{cancer_type}.m_values.csv',
                         analysis_name=analysis_name, cancer_type=cancer_types),
        beta_values=expand('results/14.{analysis_name}.{cancer_type}.beta_values.csv',
                           analysis_name=analysis_name, cancer_type=cancer_types)
    threads: nthread
    shell:
        """
        Rscript scripts/14_extractMandBetaValues.R \
            -g '{input.grsets}' \
            -b '{output.beta_values}' \
            -m '{output.m_values}'
        """


# --- 15. Collapse Adjacent CpG sites into methyalted regions for M and Beta values
rule collapse_adjacent_cpg_sites_to_methylated_regions:
    input:
        grsets=expand('results/13.{analysis_name}.{cancer_type}.GenomicRatioSet.qs', 
                        analysis_name=analysis_name, cancer_type=cancer_types)
    output:
        m_region_grsets=expand('results/15.{analysis_name}.{cancer_type}.GenomicRatioSet.m_values.regions.qs',
                                 analysis_name=analysis_name, cancer_type=cancer_types),
        beta_region_grsets=expand('results/15.{analysis_name}.{cancer_type}.GenomicRatioSet.beta_values.regions.qs',
                                  analysis_name=analysis_name, cancer_type=cancer_types)
    threads: nthread
    shell:
        """
        Rscript scripts/15_collapseAdjacentCpGSitesToMethylatedRegions.R \
            -g '{input.grsets}' \
            -b '{output.beta_region_grsets}' \
            -m '{output.m_region_grsets}'
        """


# ---- 16. Build methylated region to CpG mappings and write region M and Beta values to disk

methylation_values = ['m_values', 'beta_values']

rule build_region_mappings_for_beta_and_m_values:
    input:
        grsets=expand('results/15.{analysis_name}.{cancer_type}.GenomicRatioSet.{methylation_values}.regions.qs',
            analysis_name=analysis_name, cancer_type=cancer_types, methylation_values=methylation_values),
        methylation_data=expand('results/15.{analysis_name}.{cancer_type}.GenomicRatioSet.{methylation_values}.regions.qs',
            analysis_name=analysis_name, cancer_type=cancer_types, methylation_values=methylation_values),
        methylation_values=methylation_values
    output:
        mappings=expand('results/16.{analysis_name}.{cancer_type}.{methylation_value}.region_to_cpg_mappings.csv', 
                        analysis_name=analysis_name, cancer_type=cancer_types, methylation_value=methylation_values),
        methyl_values=expand('results/16.{analysis_name}.{cancer_type}.{methylation_value}.regions.csv', 
                             analysis_name=analysis_name, cancer_type=cancer_types, methylation_value=methylation_values)
    threads: nthread
    shell:
        """
        Rscript scripts/16_buildRegionMappingsForBetaAndMValues.R \
            -g {input.grsets} \
            -d {input.methylation_data} \
            -v {input.methylation_values} \
            -M {output.mappings} \
            -V {output.methyl_values}
        """