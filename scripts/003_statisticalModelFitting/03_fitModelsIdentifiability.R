####
####
# Study parameter identifiability in non-parametric models
####
####

# - Libraries ----
library(MCMCvis)

# - Directories ----
primaryDirectory <- "~/Dropbox/chapter-4/analysis/"
scriptDirectory <- paste0(primaryDirectory,"scripts/003_statisticalModelFitting/")
outDataDirectory <- paste0(primaryDirectory,"outputs/001_simulateObservations/01_identifiability/")
outPosteriorSamplesDirectory <- paste0(primaryDirectory,"outputs/002_statisticalModelFitting/01_identifiability/")

# - Functions ----
# - +Generate table of parameters from simulation filenames ----
f = function(x){
  tmp=strsplit(x,"-")[[1]]
  n.bags=as.numeric(tmp[2])
  p.m=as.numeric(tmp[3])
  p.g=as.numeric(sub(".RDS","",tmp[4]))
  obj=c(n.bags=n.bags,p.m=p.m,p.g=p.g)
  return(obj)
}

# - Parameters for identifiability analysis ----
# - +Build parameter table ----

fileNames <- list.files(paste0(outDataDirectory))
n  <-  length(fileNames)
parameterTable <- data.frame(do.call(rbind,lapply(fileNames,f)))

# - +Subset parameters ----
parametersIdentifiability <- parameterTable$p.m==.1 & parameterTable$p.g==0.1 & parameterTable$n.bags %in% c(100);
index = (1:n)[parametersIdentifiability]

# - Create directory ----
for(i in 1:3){

  dir.create(paste0(outPosteriorSamplesDirectory,"/posteriors-NpCmCg-",i))
  dir.create(paste0(outPosteriorSamplesDirectory,"/posteriors-NpVmCg-",i))
  
  tempDirectory <- paste0(tempDirectory,"/posteriors-NpCmCg-i")
  dir.create(paste0(tempDirectory,"/seedBagBurial"))
  dir.create(paste0(tempDirectory,"/seedAddition"))
  
  tempDirectory <- paste0(tempDirectory,"/posteriors-NpVmCg-i")
  dir.create(paste0(tempDirectory,"/seedBagBurial"))
  dir.create(paste0(tempDirectory,"/seedAddition"))
}
# 

# - Run models ----

identifiabilityBinary = TRUE

# - +Non-parametric, constant mortality
for(j in 1:length(index)){
  tmpPosteriorDirectory = outPosteriorSamplesDirectory
  
  for(t.identifiability in 1:3){
    simulation.data  <-  readRDS(paste0(paste0(outDataDirectory),fileNames[index[j]]))
    dataFull <- simulation.data[[2]]
    parms <- simulation.data[[1]]
    t.max <- max(parms$t)
    subset <- dataFull$t_obs<=t.identifiability
    data$y.s_obs <- dataFull$y.s_obs[subset]
    data$y.g_obs <- dataFull$y.g_obs[subset]
    data$n.g_obs <- dataFull$n.g_obs[subset]
    data$n.s_obs <- dataFull$n.s_obs[subset]
    data$t_obs <- dataFull$t_obs[subset]
    data$n_s <- length(data$n.s_obs)
    data$n_g <- length(data$n.g_obs)
    
    source(paste0(scriptDirectory,"01_fitNpCmCg.R"))
  }
    
}

# - +Non-parametric, varying mortality
for(j in 1:length(index)){
  
  tmpPosteriorDirectory = outPosteriorSamplesDirectory
  
  for(t.identifiability in 1:3){
    simulation.data  <-  readRDS(paste0(paste0(outDataDirectory),fileNames[index[j]]))
    dataFull <- simulation.data[[2]]
    parms <- simulation.data[[1]]
    t.max <- max(parms$t)
    subset <- dataFull$t_obs<=t.identifiability
    data$y.s_obs <- dataFull$y.s_obs[subset]
    data$y.g_obs <- dataFull$y.g_obs[subset]
    data$n.g_obs <- dataFull$n.g_obs[subset]
    data$n.s_obs <- dataFull$n.s_obs[subset]
    data$t_obs <- dataFull$t_obs[subset]
    data$n_s <- length(data$n.s_obs)
    data$n_g <- length(data$n.g_obs)
    
    source(paste0(scriptDirectory,"02_fitNpVmCg.R"))
  }
  
}
