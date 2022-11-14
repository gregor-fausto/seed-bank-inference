####
####
# Primary script to replicate analysis
####
####

# - Required packages ----

# tidyverse
# binom
# MCMCvis
# rjags

# - Create output directories ----

setwd("../")

primaryDirectory = ""

dir.create(paste0(primaryDirectory,"outputs"))
dir.create(paste0(primaryDirectory,"outputs/001_simulateObservations"))
tmp <- paste0(primaryDirectory,"outputs/001_simulateObservations")

dir.create(paste0(tmp,"/01_identifiability"))
dir.create(paste0(tmp,"/02_estimability"))
dir.create(paste0(tmp,"/03_misspecification"))

dir.create(paste0(primaryDirectory,"outputs/002_statisticalModelFitting"))
tmp <- paste0(primaryDirectory,"outputs/002_statisticalModelFitting")

dir.create(paste0(tmp,"/01_identifiability"))
dir.create(paste0(tmp,"/02_estimability"))
dir.create(paste0(tmp,"/03_misspecification"))
dir.create(paste0(tmp,"/04_overspecification"))

dir.create(paste0(primaryDirectory,"products"))
dir.create(paste0(primaryDirectory,"products/figures"))

# - Simulate observations ----

source("scripts/002_simulateObservations/01_simulateSeedBurialIdentifiability.R")
source("scripts/002_simulateObservations/03_simulateSeedBurialEstimability.R")
source("scripts/002_simulateObservations/05_simulateSeedBurialMisspecification.R")

# - Fit models to simulated observations ----

source("scripts/003_statisticalModelFitting/03_fitModelsIdentifiability.R")
source("scripts/003_statisticalModelFitting/04_fitModelsEstimability.R")
source("scripts/003_statisticalModelFitting/05_fitModelsMisspecification.R")
source("scripts/003_statisticalModelFitting/06_fitModelsOverspecification.R")

# - Analyze consequences of identifiability ----

source("scripts/004_identifiability/01_profileLogLikelihoodsConstantMortality.R")
source("scripts/004_identifiability/02_profileLogLikelihoodsAgeDependentMortality.R")
source("scripts/004_identifiability/03_identifiabilityConstantMortality.R")
source("scripts/004_identifiability/04_identifiabilityAgeDependentMortality.R")

# - Analyze statistical properties of estimates for identifiable models ----

source("scripts/005_estimability/02_analyzeEstimability.R")

# - Analyze consequences of misspecified model (C/C model fit to A/C data) ----

source("scripts/006_misspecification/02_analyzeMisspecificationMultipleGerminationValues.R")

# - Analyze consequences of misspecified model (A/C model fit to C/C data) ----

source("scripts/007_overspecification/02_analyzeOverspecification.R")

# - Run analysis for literature synthesis ----

source("scripts/008_literatureSynthesis/01_analysis.R")

# - Produce additional figures and diagrams ----

source("scripts/009_figuresDiagrams/01_conceptualFigure.R")
source("scripts/009_figuresDiagrams/02_constantMortalityConstantGermination.R")
source("scripts/009_figuresDiagrams/03_ageDependentMortalityConstantGermination.R")
