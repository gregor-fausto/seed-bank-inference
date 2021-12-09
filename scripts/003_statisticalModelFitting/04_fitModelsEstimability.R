####
####
# Study statistical properties of constant mortality, constant germination model
####
####

# - Libraries ----
library(MCMCvis)

# - Directories ----
primaryDirectory <- "../"
scriptDirectory <- paste0(primaryDirectory,"scripts/003_statisticalModelFitting/")
outDataDirectory <- paste0(primaryDirectory,"outputs/001_simulateObservations/02_estimability")
outPosteriorSamplesDirectory <- paste0(primaryDirectory,"outputs/002_statisticalModelFitting/02_estimability/")

# - Functions ----
# - +Generate table of parameters from simulation filenames ----
f = function(x){
  tmp=strsplit(x,"-")[[1]]
  n.bags=as.numeric(tmp[3])
  p.m=as.numeric(tmp[4])
  p.g=as.numeric(sub(".RDS","",tmp[5]))
  obj=c(n.bags=n.bags,p.m=p.m,p.g=p.g)
  return(obj)
}

# - Parameters for estimability ----
# - +Build parameter table ----

fileNames <- list.dirs(paste0(outDataDirectory),recursive=FALSE)
n  <-  length(fileNames)
parameterTable <- data.frame(do.call(rbind,lapply(fileNames,f)))

# filter parameter table
index = (1:n)

# create directories to hold posteriors
for(i in 1:length(index)){
  dir.create(file.path(paste0(outPosteriorSamplesDirectory,"data-",parameterTable$n.bags[i],"-",
                              parameterTable$p.m[i],"-",parameterTable$p.g[i])))
    tempDirectory <- paste0(outPosteriorSamplesDirectory,"data-",parameterTable$n.bags[i],"-",
                            parameterTable$p.m[i],"-",parameterTable$p.g[i])
    dir.create(paste0(tempDirectory,"/posteriors-NpCmCg"))
    tempDirectory <- paste0(tempDirectory,"/posteriors-NpCmCg")
    dir.create(paste0(tempDirectory,"/seedBagBurial"))
    dir.create(paste0(tempDirectory,"/seedAddition"))
}

# - Run models ----

identifiabilityBinary = FALSE

# - +number of replicate simulations ----
n.replicate = 10

simulatedData = paste0("replicate-",1:n.replicate,".RDS")


# fit non-parametric seed bag burial and seed addition experiments
for(i in 1:length(index)){

  simulatedDataObj <- strsplit(fileNames[index[i]],"/")[[1]][10]

  tmpPosteriorDirectory <- paste0(outPosteriorSamplesDirectory,simulatedDataObj,"/")

  for(j in 1:n.replicate){

    # if replicate exists go on
    if(file.exists(paste0(tmpPosteriorDirectory,"posteriors-NpCmCg/seedBagBurial/replicate-",j,".RDS"))) next
    simulation.data  <-  readRDS(paste0(paste0(fileNames[index[i]],"/"),simulatedData[j]))
    source(paste0(scriptDirectory,"01_fitNpCmCg.R"))

  }
}
