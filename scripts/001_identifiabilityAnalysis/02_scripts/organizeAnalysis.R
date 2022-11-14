####
####
# Compile and organize results of identifiability analysis
# these results will be the basis for Table 2
####
####

# - Libraries ----
library(tidyverse)

# - Read data ----
# - +file names ----
fileNames <- read.table("scripts/001_identifiabilityAnalysis/03_output/filenames.txt",sep=";")
# - +model deficiency ----
deficiency <- read.table("scripts/001_identifiabilityAnalysis/03_output/symbolicOutput.txt",sep=";")

# - Extract variables from file names ----

# - +read filename as string and split by backslash ----
fileNames <- as.character(fileNames[,1])
fileNamesSplit<-strsplit(fileNames,"/")

# - +experiment ----
# note that the goal here is to extract either 'seed addition' or 'seed bag burial' 
# may need to change the index from 10 to whatever matches that variable
# if the file name changes
experiment <-unlist(lapply(fileNamesSplit, `[[`, 10))

# - +model structure ----
# may need to change the index from 11 
# if the file name changes
# you want to recover the variable that has for eg. 'gc-ma-1'
model <-unlist(lapply(fileNamesSplit, `[[`, 11))

# - Extract model deficiency ----
deficiency <- deficiency[,1]

# - Bind experiment, model, and deficiency in data frame ----
deficiencyResults <- data.frame(experiment = experiment, model = model, deficiency = deficiency) 

# - +Reorganize data frame ----
deficiencyResults <- deficiencyResults %>%
  dplyr::mutate(model=as.character(model)) %>%
  dplyr::mutate(years = substring(model,7)) %>%
  dplyr::mutate(model = substring(model,1,5))

# - +Reshape data frame for table in paper ----
deficiencyResults <- deficiencyResults %>%
  tidyr::spread("years","deficiency",deficiency)

# - +Create data frame with experiment and variable names ----
analysisVariables <- data.frame(experiment=c(rep("seedAddition",4),rep("seedBagBurial",4)),
                     model = c("gc-mc","gc-ma","ga-mc","ga-ma",
                               "gc-mc","gc-ma","ga-mc","ga-ma"))

# - +join data frames ----
identifiabilityOutput <- left_join(analysisVariables,deficiencyResults,by=c("experiment","model"))

# - Print output ----
identifiabilityOutput
