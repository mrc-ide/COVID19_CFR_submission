# COVID19_CFR_submission
Repository for all scripts required to replicate the analysis in the paper ["Estimates of the severity of coronavirus disease 2019: a model-based analysis"](https://doi.org/10.1016/S1473-3099(20)30243-7).

All code within this repository is freely available and published under the MIT license (see LICENSE file)

## How to run the code

Code in this repository is written in the [R statistical language](https://www.r-project.org/), and makes use of [RStudio projects](https://rstudio.com/) for convenience along with several additional packages. For those just starting out in R check out the [swirl](https://swirlstats.com/) package, and for more information on how to work with packages check out [this page](https://www.datacamp.com/community/tutorials/r-packages-guide).

We recommend starting by opening "COVID19_CFR_submission.Rproj" within RStudio, this will ensure the working directory is set correctly. The various scripts required to run the analyses are contained in the R_scripts folder. These scripts were written by multiple authors, and so style and implementation differs between scripts, but in general these scripts will install/load any required packages and functions at the top of the script, the body of the script then carries out the main analysis and any output is saved to the "output" folder.

Scripts required to replicate the parametric analysis in Table 2 begin with "run_mcmc_cfr", and can be run in any order. These scripts run MCMC via the [drjacoby package](https://github.com/mrc-ide/drjacoby). They make use of processed data in the "output" folder. If you want to see the raw data prior to this processing step (this is not required for the analysis) check out the "data" folder and the script "process_mcmc_data.R".


