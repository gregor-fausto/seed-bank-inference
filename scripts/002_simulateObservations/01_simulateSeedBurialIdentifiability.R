####
####
# Simulate observations for a range of parameter values
####
####

# - Set seed ----
set.seed(16)

# - Directories ----
primaryDirectory <- ""
scriptDirectory <- paste0(primaryDirectory,"scripts/002_simulateObservations/")
outDataDirectory <- paste0(primaryDirectory,"outputs/001_simulateObservations/01_identifiability/")

# - Parameters for simulation ----

# - +Probability of germination ----
p.g <- 0.1

# - +Probability of mortality ----
p.m <- 0.1

# - +Number of seeds per bag/plot at start of experiment ----
n.start <- 100

# - +Number of bags/plots ----
n.bags <- 1000

# - Build table ----
parameterTable <- expand.grid(p.g=p.g,p.m=p.m,n.start=n.start,n.bags=n.bags)

# filter to unique combinations
parameterTable <- unique(parameterTable)

# - Simulate observations ----
n=dim(parameterTable)[1]

# - +for each row of the parameter table ----
for(j in 1:n){
  tmp <- parameterTable[j,]

  # - +extract parameters ----
  p.g <- tmp$p.g;p.m <- tmp$p.m;n.start <- tmp$n.start;n.bags <- tmp$n.bags
  # - +set time steps to 1, 2, 3 years ----
  t <- 1:3
  # - +source script for simulation ----
  source(paste0(scriptDirectory,"00_functionSimulateSeedBurialIdentifiability.R"))
  
  # - +save output ----
  saveRDS(list(parameters,data),file=paste0(paste0(outDataDirectory),"data-",n.bags,"-",p.m,"-",p.g,".RDS"))
}
