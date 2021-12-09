####
####
# Primary script to replicate analysis
####
####

# - Simulate observations ----

source("002_simulateObservations/01_simulateSeedBurialIdentifiability.R")
source("002_simulateObservations/03_simulateSeedBurialEstimability.R")
source("002_simulateObservations/05_simulateSeedBurialMisspecification.R")

# - Fit models to simulated observations ----

source("003_statisticalModelFitting/03_fitModelsIdentifiability.R")
source("003_statisticalModelFitting/04_fitModelsEstimability.R")
source("003_statisticalModelFitting/05_fitModelsMisspecification.R")

# - Analyze consequences of identifiability ----

source("004_identifiability/01_profileLogLikelihoodsConstantMortality.R")
source("004_identifiability/02_profileLogLikelihoodsAgeDependentMortality.R")
source("004_identifiability/03_identifiabilityConstantMortality.R")
source("004_identifiability/04_identifiabilityAgeDependentMortality.R")

# - Analyze statistical properties of estimates for identifiable models ----

source("005_estimability/02_analyzeEstimability.R")

# - Analyze consequences of misspecified model ----

source("006_misspecification/01_analyzeMisspecification.R")

# - Run analysis for literature synthesis ----

source("007_literatureSynthesis/01_analysis.R")

# - Produce additional figures and diagrams ----

source("008_figuresDiagrams/01_conceptualFigure.R")
source("008_figuresDiagrams/02_constantMortalityConstantGermination.R")
source("008_figuresDiagrams/03_ageDependentMortalityConstantGermination.R")
