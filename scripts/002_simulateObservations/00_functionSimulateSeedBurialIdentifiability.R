####
####
# Simulate observations from an experiment
####
####

# - Function to obtain hazards based on annual probability of germination and mortality ----
f.hazards = function(t= t, p.g=.1,p.m=.1){

  # get discrete time point
  time.discrete = unique(floor(t))

  ## Set germination hazards
  # If 1 value is supplied for p.g, assume single hazard at all times
  # If fewer values supplied for p.g than time points, extend last hazard to end of time
  # If number of values for p.g matches time points, set p.g equal to hazard
  n.t = length(time.discrete)#-1
  if(length(p.g) < n.t){
    if(length(p.g)==1){
      hazard.g = c()
      hazard.g=rep(p.g,n.t)
    } else {
      n.tmp=length(p.g)
      hazard.g=c(p.g,rep(p.g[n.tmp],n.t-n.tmp))
    }
  } else {
    hazard.g = p.g
  }

  ## Set mortality hazards
  # If 1 value is supplied for p.m, assume single hazard at all times
  # If fewer values supplied for p.m than time points, extend last hazard to end of time
  # If number of values for p.m matches time points, set p.g equal to hazard
  n.t = length(time.discrete)#-1
  if(length(p.m) < n.t){
    if(length(p.m)==1){
      hazard.m = c()
      hazard.m=rep(p.m,n.t)
    } else {
      n.tmp=length(p.m)
      hazard.m=c(p.m,rep(p.m[n.tmp],n.t-n.tmp))
    }
  } else {
    hazard.m = p.m
  }

  return(list(hazard.m,hazard.g))
}

# - Function to obtain survival function based on annual probability of germination and mortality ----
f.survival = function(t=t,h.g=.2,h.m=.2){

  # get discrete time point
  time.discrete = unique(floor(t))

  # mortality component
  s.m = (1-h.m)
  s.m = c(s.m)
  # if fractional, less than integer
  if(t%%1!=0){index = time.discrete+1} else { index = time.discrete}
  s.mortality = cumprod(s.m)[1:(time.discrete+1)][index]

  # germination component
  s.g = (1-h.g)
  s.g = c(1,s.g)
  s.germination = cumprod(s.g[1:(time.discrete+1)])[time.discrete+1]

  survival.function = s.mortality*s.germination

  return(survival.function)
}

# -  Function to obtain probability mass function for germination ----
f.pmf = function(t=time.discrete,h.g=.1,h.m=.1){

  # Probability mass function for emergence/recruitment
  surv=sapply(t-exp(-26),f.survival,h.g,h.m)
  p.emergence = c()
  for(i in t){
    p.emergence[i] = surv[i]*h.g[i]
  }

  p.mortality = c()
  n = length(time.discrete)
  int.start=c(1,(surv-p.emergence)[1:(n-1)])
  p.mortality = int.start-surv

  return(list(p.mortality,p.emergence))
}


# - Simulation for a seed bank in which mortality precedes germination ----
# Discrete times

# - +Calculate hazards ----
h.m=f.hazards(t=t,p.g=p.g,p.m=p.m)[[1]]
h.g=f.hazards(t=t,p.g=p.g,p.m=p.m)[[2]]

survival.yearend = sapply(t-exp(-16),f.survival,h.g=h.g,h.m=h.m)

# - +Calculate probability mass functions ----
time.discrete=unique(floor(t))

pmf.m = f.pmf(t=time.discrete,h.g=h.g,h.m=h.m)[[1]]
pmf.g = f.pmf(t=time.discrete,h.g=h.g,h.m=h.m)[[2]]

# - +Simulate observations ----
y.g.sample=matrix(NA,nrow=n.bags,ncol=length(t))
y.s.sample=matrix(NA,nrow=n.bags,ncol=length(t))

for(i in t){
  p.s.tmp=survival.yearend
  y.s.sample[,i]=rbinom(n.bags,n.start,prob=p.s.tmp[i])
  n.g.sample=y.s.sample
  y.g.sample[,i]=rbinom(n.bags,n.g.sample[,i],prob=h.g[i])
}

# - +Build list ----
y.s_obs = as.vector(y.s.sample)
n.s_obs = rep(n.start,length(y.s_obs))
t_obs = rep(t,each=n.bags)
y.g_obs = as.vector(y.g.sample)
n.g_obs = as.vector(n.g.sample)
n_s = length(y.s_obs)
n_g = length(y.g_obs)

data = list(y.s_obs=y.s_obs,y.g_obs=y.g_obs,
            n.g_obs=n.g_obs,n.s_obs=n.s_obs,
            n_s=n_s,n_g=n_g,t_obs=t_obs)

parameters = list(t=t,p.g=p.g,p.m=p.m,n.start=n.start,n.bags=n.bags)

