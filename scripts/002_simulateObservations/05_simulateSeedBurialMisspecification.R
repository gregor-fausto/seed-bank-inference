####
####
# Simulate observations from an age-dependent mortality, constant germination process
####
####

# - Set seed ----
set.seed(16)

# - Directories ----
primaryDirectory <- "~/Dropbox/chapter-4/analysis/"
scriptDirectory <- paste0(primaryDirectory,"scripts/002_simulateObservations/")
outDataDirectory <- paste0(primaryDirectory,"outputs/001_simulateObservations/03_misspecification/")

# - Parameters for simulation ----

# - +Probability of germination----
p.g <- c(.1)

# - +Probability of mortality in first year----
p.m <- c(.1)

# - +Number of seeds per bag/plot at start of experiment ----
n.start <- 100

# - +Number of bags/plots ----
n.bags <- seq(5,30,by=5)

parameterTable <- expand.grid(p.g=p.g,p.m=p.m,n.start=n.start,n.bags=n.bags)

# - Build table ----
parameterTable=rbind(parameterTable)

# - Simulate observations ----
n=dim(parameterTable)[1]

# - + Set number of replicate simulated datasets ----
n.replicate = 10

# - +for each row of the parameter table ----
for(h in 1:n){
  
  tmp <- parameterTable[h,]
  
  # - +extract parameters ----
  p.g <- tmp$p.g;p.m <- tmp$p.m;n.start <- tmp$n.start;n.bags <- tmp$n.bags
  
  # - +modify mortality so that it increases by 0.1 each year ----
  p.m = c(p.m,p.m+.1,p.m+.2)
  
  # - +set time steps to 1, 2, 3 years ----
  t <- 1:3
  
  # - +for each parameter set, create n.replicate simulated datasets ----
  for(j in 1:n.replicate){
    
    # - +source script for simulation ----
    source(paste0(scriptDirectory,"04_functionSimulateSeedBurialMisspecification.R"))
    
  }
}

