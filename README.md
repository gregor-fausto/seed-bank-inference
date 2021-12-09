# seed-bank-inference
Scripts associated with paper on seed bank inference

### Statistical inference for seed mortality and germination with seed bank experiments

### Authors

  - Gregor-Fausto Siegmund, Cornell University, <gs589@cornell.edu>
  - Monica A. Geber, Cornell University

  This repository contains the analysis for our project on statistical models for seed mortality and germination.

-----

### Abstract

….

-----

### Repository Directory

The repository is organized so that the analysis in the paper can be replicated.

- `models`: Contains the statistical models written in JAGS language.
    + `seedAddition`: models for seed addition experiments
    + `seedBagBurial`: models for seed bag burial experiments
- `scripts`: Contains R scripts to run the simulations, fit models, and analyze output.
    + `001_identifiabilityAnalysis`: scripts and functions for algebraic identifiability analysis
    + `002_simulateObservations`: scripts to simulate observations from seed bag burial and seed addition experiments
    + `003_statisticalModelFitting`: scripts to fit JAGS models to simulated observations
    + `004_identifiability`: scripts to analyze effect of identifiability on parameter estimation
    + `005_estimability`: scripts to compare statistical properties of estimates from seed bag burial vs. seed addition experiments
    + `006_misspecification`: scripts to compare estimates from model with incorrect assumptions for seed bag burial vs. seed addition experiments
    + `007_literatureSynthesis`: script to analyze literature synthesis
    + `008_figuresDiagrams`: scripts to create additional figures and diagrams for manuscript
- `outputs`: Folder for output of simulations and model fitting
    + `001_simulateObservations`: simulated observations from seed bag burial  and seed addition experiments
        * `01_identifiability`: simulated observations for showing effect of identifiability on parameter estimation.         
        * `02_estimability`: simulated observations to compare statistical properties of estimates from seed bag burial vs. seed addition experiments
        * `03_misspecification`: simulated observations to evaluate consequences of incorrect assumptions about seed mortality
    + `002_statisticalModelFitting`: scripts to fit JAGS models to simulated observations
        * `01_identifiability`: samples from posterior distribution of model fits         
        * `02_estimability`: samples from posterior distribution of model fits
        * `03_misspecification`: samples from posterior distribution of model fits
- `data`: Contains data collected for the study
    + `001_literatureSynthesis`: CSV file with data collected for the literature synthesis
- `products`:
