# -------------------------------------------------------------------
# Simulation
# -------------------------------------------------------------------

data <- simulation.data[[2]]
parms <- simulation.data[[1]]
t.max <- max(parms$t)

# -------------------------------------------------------------------
# Set JAGS parameters and random seed
# -------------------------------------------------------------------
# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain

n.chain <- 3
n.adapt <- 3000
n.update <- 5000
n.iterations <- 5000
n.thin <- 10

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Negative exponential, constant germination
# -------------------------------------------------------------------
# -------------------------------------------------------------------
library(rjags)

# set inits for JAGS
inits.fun <- function(n.pars = 1){
  rbeta(n = n.pars, 1, 1 )
}


inits <- list()

for(init.i in 1:3){
  inits[[init.i]] <- list(inits.fun(n.pars=3) , inits.fun(n.pars=1))
  
  names(inits[[init.i]]) = c("p.m","p.g")
  
}

# tuning (n.adapt)
jm <- jags.model(paste0(primaryDirectory,"models/jags/seedAddition/jagsNpVariableMortalityConstantGermination.R"), data = data, n.adapt = n.adapt,
                 n.chains=3)

# burn-in (n.update)
update(jm, n.iter = n.update)

samples.rjags <- coda.samples(jm,
                              variable.names = c("p.m","p.g"),
                              n.iter = n.iterations, thin = n.thin)

if(identifiabilityBinary == TRUE){ 
  saveRDS(samples.rjags,file=paste0(outPosteriorSamplesDirectory,"posteriors-NpVmCg-",t.identifiability,"/seedAddition/posterior-",parms$n.bags,"-",parms$p.m,"-",parms$p.g,".RDS"))
} else if(identifiabilityBinary == FALSE){
  saveRDS(samples.rjags,file=paste0(tmpPosteriorDirectory,"posteriors-NpVmCg/seedAddition/replicate-",j,".RDS"))
}

rm(list=setdiff(ls(all=TRUE),
                c("primaryDirectory","data","n.adapt","n.update", "n.chain", "n.iterations", "n.thin","j","i","tmpPosteriorDirectory","n.sims",
                  "outPosteriorSamplesDirectory", "outDataDirectory", "parms", "t.max","fileNames","index","scriptDirectory","simulatedData",
                  "t.identifiability","identifiabilityBinary","n.replicate")))

library(rjags)

# set inits for JAGS
inits.fun <- function(n.pars = 1){
  rbeta(n = n.pars, 1, 1 )
}


inits <- list()

for(init.i in 1:3){
  inits[[init.i]] <- list(inits.fun(n.pars=3) , inits.fun(n.pars=1))
  
  names(inits[[init.i]]) = c("p.m","p.g")
  
}


# tuning (n.adapt)
jm <- jags.model(paste0(primaryDirectory,"models/jags/seedBagBurial/jagsNpVariableMortalityConstantGermination.R"), data = data, n.adapt = n.adapt,
                 n.chains=3)



# burn-in (n.update)
update(jm, n.iter = n.update)

samples.rjags <- coda.samples(jm,
                              variable.names = c("p.m","p.g"),
                              n.iter = n.iterations, thin = n.thin)

if(identifiabilityBinary == TRUE){ 
  saveRDS(samples.rjags,file=paste0(outPosteriorSamplesDirectory,"posteriors-NpVmCg-",t.identifiability,"/seedBagBurial/posterior-",parms$n.bags,"-",parms$p.m,"-",parms$p.g,".RDS"))
} else if(identifiabilityBinary == FALSE){
  saveRDS(samples.rjags,file=paste0(tmpPosteriorDirectory,"posteriors-NpVmCg/seedBagBurial/replicate-",j,".RDS"))
}

rm(list=setdiff(ls(all=TRUE),
                c("primaryDirectory","data","n.adapt","n.update", "n.chain", "n.iterations", "n.thin","j","i","tmpPosteriorDirectory","n.sims",
                  "outPosteriorSamplesDirectory", "outDataDirectory", "parms", "t.max","fileNames","index","scriptDirectory","simulatedData",
                  "t.identifiability","identifiabilityBinary","n.replicate")))
