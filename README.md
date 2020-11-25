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
the R scripts.

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

Please at minimum create the `rawdata` and `metadata` directories, as the
are assumed to hold the raw Illumina EPIC microarray plate data (.xml, .idat, etc.) and the plate labels, respectively. The remaining missing directories will be
created automatically as the pipeline runs.

### config.yaml

This file hold the relevant pipeline documentation.