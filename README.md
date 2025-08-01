# The Winners Take It All? Evolutionary Success of H5Nx Reassortants in the 2020–2024 Panzootic

This repository contains code used for the data handling and statistical 
analysis of Clade 2.3.4.4b High Pathogenicity Avian Influenza Virus (HPAIV) 
circulating during the 2020-2024 Panzootic. Specifically, we analysed the emergence,
persistence and drivers of unique reassortants worldwide, up to May 2024. 

## Data Curation



## Statistical Models
We fitted three statistical models to quantify patterns of reassortant emergence
across continents and to understand the drivers of reassortant spatial diffusion. 
All models are located within the [statistical_models](scripts/statistical_models/) 
sub directory.\\

Briefly, the 'number of reassortants model' is a mixture model comprising a 
zero-inflated Poisson process and a binomial 'filter' process; the 'reassortant class model'
is an ordinal model, assuming a cumulative distribution; and the 'diffusion model' 
is a log-gamma regression. The
'number of reassortants model' was fitted using Stan via rStan, and the remainder were
fitted using BRMS. Each model has three associated scripts: one 
fitting the model, one to run model evaluations and one for interpretation. The
'number of reassortants model' has additional scripts to describe the model in Stan 
(located in [scripts/statistical_models/stan_models](scripts/statistical_models/stan_models))
and for pre-processing.

Evaluation and interpretation plots are produced in \*_model_evaluation and \*model_interpretation
scripts, however these may differ from the final published plots (located in
[scripts/figure_scripts](scripts/figure_scripts))

All models have been tested on i) Apple M4 Max 16-core CPU with 48 Gb RAM 
and AMD Ryzen 9 7950X 16 Core CPU with 96 Gb RAM.