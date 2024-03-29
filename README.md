# seed-bank-inference

### Statistical inference for seed mortality and germination with seed bank experiments

### Authors

  - Gregor-Fausto Siegmund, Cornell University, <gs589@cornell.edu>
  - Monica A. Geber, Cornell University

  This repository contains the scripts for our project on statistical models for seed mortality and germination.

  [![DOI](https://zenodo.org/badge/436662955.svg)](https://zenodo.org/badge/latestdoi/436662955)

-----

### Abstract

Plant population ecologists regularly study soil seed banks with seed bag burial and seed addition experiments. These experiments contribute crucial data to demographic models, but we lack standard methods to analyze them. Here, we propose statistical models to estimate seed mortality and germination with observations from these experiments. We develop these models following principles of event history analysis, and analyze their identifiability and statistical properties by algebraic methods and simulation. We demonstrate that seed bag burial, but not seed addition experiments, can be used to make inferences about age-dependent mortality and germination. When mortality and germination do not change with seed age, both experiments produce unbiased estimates but seed bag burial experiments are more precise. However, seed mortality and germination estimates may be inaccurate when the statistical model that is fit makes incorrect assumptions about the age-dependence of mortality and germination. The statistical models and simulations that we present can be adopted and modified by plant population ecologists to strengthen inferences about seed mortality and germination in the soil seed bank.

-----

### Repository Directory

The repository is organized so that the analysis in the paper can be replicated.

- `models`: Contains the statistical models written in JAGS language.
    + `seedAddition`: models for seed addition experiments
    + `seedBagBurial`: models for seed bag burial experiments
- `data`: Contains data collected for the study
    + `001_literatureSynthesis`: CSV file with data collected for the literature synthesis
- `scripts`: Contains R scripts to run the simulations, fit models, and analyze output.
    + `001_identifiabilityAnalysis`: scripts and functions for algebraic identifiability analysis
    + `002_simulateObservations`: scripts to simulate observations from seed bag burial and seed addition experiments
    + `003_statisticalModelFitting`: scripts to fit JAGS models to simulated observations
    + `004_identifiability`: scripts to analyze effect of identifiability on parameter estimation
    + `005_estimability`: scripts to compare statistical properties of estimates from seed bag burial vs. seed addition experiments
    + `006_misspecification`: scripts to compare estimates from model with incorrect assumptions for seed bag burial vs. seed addition experiments; C/C model fit to A/C data
    + `007_overspecification`: scripts to compare estimates from model with incorrect assumptions for seed bag burial vs. seed addition experiments; A/C model fit to C/C data
    + `008_literatureSynthesis`: script to analyze literature synthesis
    + `009_figuresDiagrams`: scripts to create additional figures and diagrams for manuscript
    + `primaryScript`: script to run the analysis and produce figures

Running `primaryScript` in the appropriate directory will create the folders `outputs` and `products` with the following file structure. Note that replicating the simulation and model fitting may be slow. We recommend testing the code in `003_statisticalModelFitting` on a smaller number of replicates than the default, which we set as 250 replicates as in the text of the associated manuscript.

- `outputs`: Folder for output of simulations and model fitting
    + `001_simulateObservations`: simulated observations from seed bag burial  and seed addition experiments
        * `01_identifiability`: simulated observations for showing effect of identifiability on parameter estimation.         
        * `02_estimability`: simulated observations from a seed bank with constant mortality/constant germination; used to compare statistical properties of estimates from seed bag burial vs. seed addition experiments and for analysis of over-specification (fitting the A/C model to C/C data)
        * `03_misspecification`: simulated observations from a seed bank with age-dependent mortality/constant germination; used to evaluate consequences of incorrect assumptions about seed mortality (fitting the C/C model to A/C data)
    + `002_statisticalModelFitting`: scripts to fit JAGS models to simulated observations
        * `01_identifiability`: samples from posterior distribution of model fits         
        * `02_estimability`: samples from posterior distribution of model fits
        * `03_misspecification`: samples from posterior distribution of model fits
        * `04_overspecification`: samples from posterior distribution of model fits
- `products`: Folder for figures
    + `figures`: directory to hold figures produced by scripts      
