import re

## TODO: Implement check to see if packrat dependencies are available and install if not

configfile: 'config.yaml'

plate_data_dirs = config['plate_data_dirs']
plate_labels = config['plate_labels']
nthread = config['nthread']
analysis_name = config['analysis_name']

# ---- 0.


# ---- 1. Build minfi RGSet from raw plate data and labels

rule build_rgset_from_plate_data:
    input:
        plates=expand('rawdata/{plate_dirs}', plate_dirs=plate_data_dirs),\
        labels=expand('metadata/{labels}', labels=plate_labels)
    output:
        expand('rawdata/1_{analysis_name}_RGSet_raw.qs', analysis_name=analysis_name)
    shell:
        """
        Rscript scripts/buildRGsetFromPlateData.R -p '{input.plates}' -l '{input.labels}' -o {output} -n {nthread}
        """

# ---- 2. 