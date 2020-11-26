# illuminaEPICmethylation

This repository contains a pipeline utilizing the `minfi` R package to 
process raw Illumina Infinum 850k array data into sample by feature .csv 
files containing M and Beta value calls. While this pipeline was 
specifically developed to analyze 850k arrays, it could easily be adapted 
for use with 450k arrays by changing some of the reference and annotation
packages loaded in R scripts used in the pipeline.

In the future, we may consider parameterizing libraries to allow this
pipeline to support the full range of potential `minfi` inputs, but for
the time being users will need to manually change the library calls in
the provided R scripts.

## Snakemake

This pipeline leverages the `snakemake` Python package for workflow management. As a result the pipeline and its dependencies are easily
installed from this repository, allowing quick setup and configuration
before the pipeline can be deployed.

For more information on Snakemake, please see: https://snakemake.readthedocs.io/en/stable/.

## Requirements

Dependency management for this pipeline is handled via `conda` for Python 
and `renv` for R. To get started with setup you can install
miniconda3 using the instructions available here: https://docs.conda.io/en/latest/miniconda.html. If you do not currently have R installed, you can install it via conda using the command: `conda -c conda-forge r-base`. 

Alternatively you can install it directly from CRAN
as described here: https://cran.r-project.org/.

## Setting Up Your Software Environment

The first step to deploying an analysis pipeline is to install the various
software packages it depends on. We have included the `envs/environment.yml` and `renv.lock` files here to easily accomplish this.

All commands should be executed from the top level directory of this
repository.

### Python and System Dependencies

Conda can be used to install various Python and system dependencies
using:

`conda env create --file envs/environment.yml`

This will take some time to run as it gathers and installs the correct
package versions. The environent it creates should be called `methylation`.

If it is not automatically activated after installation please run: 
`conda activate methylation` before proceeding to the next step.

### R Dependencies

The `renv` package can be used to install all R dependencies (both CRAN and
Bioconductor). Please ensure you have R version >=4.0 installed. The 
pipeline may work with older versions or R, but it has not been tested.

To initialize this project with renv, first open R interactively and run:
`install.packages('renv')`. This is required since you need to select a
CRAN mirror to install from.

Once installed, you can exit R and run:

`Rscript -e 'library(renv); renv::init()'`

If you wish to isolate the R dependencies from your system wide R libraries, you can use this command instead:

`Rscript -e 'library(renv); renv::isolate(); renv::init()'`

If intialization doesn't trigger dependency installation, you can do so manually using:

`Rscript -e 'renv::restore()'`

For more information on renv and it can be used to manage dependencies in
your project, please see: https://rstudio.github.io/renv/articles/renv.html.

## Configuring the Pipeline

This pipeline assumes the following directory structure:

```
.
├── envs
├── metadata
├── procdata
├── qc
├── rawdata
├── renv
├── results
└── scripts
```

Please at minimum create the `rawdata` and `metadata` directories, as they
are assumed to hold the raw Illumina EPIC microarray plate data (.xml, .idat, etc.) and the plate labels, respectively. The remaining missing directories will be
created automatically as the pipeline runs.

### config.yaml

This file hold the relevant pipeline documentation. At minimum, you
must specify the paths to the raw plate data to `plate_data_dirs`,
the plate labels to `plate_label`, the max `nthread` the pipeline should
use, and an `analysis_name`. See the comments in the file for more information.

## Using the Pipeline

### Read in Data and Do Plate Level QC

If you open `Snakefile` you will see the details of the workflow that this pipeline
executes. The steps of the pipeline are organizeed into rule, and the output file
from each rule is prepended with the rule number and the `analysis_name`.

Because there are several qualitative QC steps in this pipeline, we do not recommend
running it end to end. Once you configure the parameters in `config.yaml` up to
`manual_qc_steps`, you should run the command:

`snakemake --cores 2 drop_samples_with_bisulphite_conversion_less_than_cutoff`

This will execute rule 1 through 3 using the number of cores specfied to `--cores`. 

Step 1 creates an RGChannelSet from the raw plate data. Step 2 filters out samples
with less the 90% of CpG probes detected, it also generates a number of qc files
in the `qc` directory. Step 3 filters out probes with a bisulphite conversion rate
lower than what is specified in `config.yaml`.

At this point you should have a look at the the QC metrics generated in the previous step
and decide which, if any, samples you wish to filter out. There may also be samples you
wish to remove for other reasons, such as techincal replicates or wells spiked for qc.

To do this you will use specify one or more rounds of sample removal to the `manual_qc_steps`
parameter in `config.yaml`. Please see the documentation provided in that file.

### Manual QC and Preprocessing

Once you have setup your manual qc steps, we recommend running the command:

`snakemake --cores 2 density_plot_preprocessed_vs_rgset`

This will execute rules 4 through 7. This rule generates two qc files. 
The first is a `.pdf` file showing the beta-value distribution for each sample under each of the normalization methods specified to the `preprocess_methods` parameter in `config.yaml`. It is a good idea to have a look at these files to determine the correct parameter cut-offs for the second round of manual qc. You can also identify samples which technically fail one of the three qc metrics computed, but still have a good distribution on visual inspection of the plot.


### Final QC and Cancer Type Separation



