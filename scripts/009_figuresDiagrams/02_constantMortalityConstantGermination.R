
# -------------------------------------------------------------------
# Function to obtain hazards based on annual probability of germination and mortality
# -------------------------------------------------------------------
f.hazards = function(t= time.discrete, p.g=.1,p.m=.1){
  
  ## Set germination hazards
  # If 1 value is supplied for p.g, assume single hazard at all times
  # If fewer values supplied for p.g than time points, extend last hazard to end of time
  # If number of values for p.g matches time points, set p.g equal to hazard
  if(length(p.g)<length(t)){
    if(length(p.g)==1){
      hazard.g = c()
      hazard.g=rep(p.g,length(t)) 
    } else {
      n.tmp=length(p.g)
      hazard.g=c(p.g,rep(p.g[n.tmp],length(t)-n.tmp))
    }
  } else {
    hazard.g = p.g
  }
  
  ## Set mortality hazards
  # If 1 value is supplied for p.m, assume single hazard at all times
  # If fewer values supplied for p.m than time points, extend last hazard to end of time
  # If number of values for p.m matches time points, set p.g equal to hazard
  if(length(p.m)<length(t)){
    if(length(p.m)==1){
      hazard.m = c()
      hazard.m=rep(p.m,length(t)) 
    } else {
      n.tmp=length(p.m)
      hazard.m=c(p.m,rep(p.m[n.tmp],length(t)-n.tmp))
    }
  } else {
    hazard.m = p.m
  }
  return(list(hazard.m,hazard.g)) 
  
}

# -------------------------------------------------------------------
# Function to obtain survival function based on annual probability of germination and mortality
# -------------------------------------------------------------------
f.survival = function(t=time.discrete,h.g=.1,h.m=.1){
  
  survival.fun = 1
  index.g=index.d=c()
  index.surv = seq(1,length(t)*2,by=2)
  for(i in t){
    tmp = index.surv[i]
    survival.fun[tmp+1] = survival.fun[tmp]*(1 - h.m[i])
    survival.fun[tmp+2] = survival.fun[tmp+1]*(1-h.g[i])
    index.g[tmp] = tmp+2;index.d[tmp]=tmp+1
  }
  
  return(survival.fun)
}

# -------------------------------------------------------------------
# Function to obtain probability mass function for germination
# -------------------------------------------------------------------
f.pmf = function(t=time.discrete,h.g=.1,h.m=.1){
  
  # Probability mass function for emergence/recruitment
  survival.fun = 1
  p.emergence = c()
  index.g=index.d=c()
  for(i in t){
    p.emergence[i] = (survival.fun[i]*(1 - h.m[i]))*h.g[i]
    survival.fun[i+1] = (survival.fun[i]*(1 - h.m[i]))*(1-h.g[i])
  }
  
  survival.fun = 1
  p.mortality = c()
  index.g=index.d=c()
  for(i in t){
    p.mortality[i] = (survival.fun[i]*(h.m[i]))
    survival.fun[i+1] = (survival.fun[i]*(1 - h.m[i]))*(1-h.g[i])
  }
  
  return(list(p.mortality,p.emergence))
}


# -------------------------------------------------------------------
# Simulation for a seed bank in which mortality precedes germination
# Discrete times
# -------------------------------------------------------------------

# define a vector for discrete times as close to integer as possible
time.discrete = 1:5

# Probability of germination
p.g = .25

# Probability of mortality
p.m = list(c(0),c(.2),c(.4))

h.m=h.g=survival.fun=pmf.m=pmf.g=list()

for(i in 1:3){
# Calculate hazards
h.m[[i]]=f.hazards(t=time.discrete,p.g=p.g,p.m=p.m[[i]])[[1]]
h.g[[i]]=f.hazards(t=time.discrete,p.g=p.g,p.m=p.m[[i]])[[2]]

# Calculate survival function
survival.fun[[i]] = f.survival(t=time.discrete,h.g=h.g[[i]],h.m=h.m[[i]])

# Calculateprobability mass functions
pmf.m[[i]] = f.pmf(t=time.discrete,h.g=h.g[[i]],h.m=h.m[[i]])[[1]]
pmf.g[[i]] = f.pmf(t=time.discrete,h.g=h.g[[i]],h.m=h.m[[i]])[[2]]
}

# -------------------------------------------------------------------
# Define function to create set of plots
# -------------------------------------------------------------------

plot.colors = c("#1b9e77","#d95f02","#7570b3")
plot.shapes = c(19,21,4)
plot.lines = c("solid","dashed","dotted")

pdf("~/Dropbox/chapter-4/analysis/products/figures/figure-2-constant.pdf",height=6,width=4)
par(mfrow=c(2,1))
par(mar=c(1,4,0,2),
    oma=c(4,0,2,0))

# Germination hazard (age-specific)

plot(NA,NA, 
     ylim=c(0,1), xlim=c(0,max(time.discrete)),
     xaxt='n',xlab="",
     axes=FALSE,
     ylab = "Germination hazard",
     pch=19)

for(i in 1:3){
  points(time.discrete,h.g[[i]],col="black",pch=19,type='b')
}

axis(1L,labels=FALSE);axis(2L)

# Mortality hazard
plot(NA,NA,
     ylim=c(0,1),xlim=c(0,max(time.discrete)),
     pch=19, axes=FALSE,
     ylab= "Mortality hazard")
for(i in 1:3){
  points(time.discrete,h.m[[i]],pch=19,col=plot.colors[i],type='b',lty=plot.lines[i])
}
axis(1L);axis(2L)
mtext("Time (t)", side = 1,line=2.2)


par(mfrow=c(1,1))
par(mar=c(1,4,0,2),
    oma=c(4,0,2,0))

# Calculations for plotting
offset = .1
time.points=rep(time.discrete,each=2) + rep(c(-offset,offset),length(time.discrete))
time.points=c(0,time.points)

km.points=rep(time.discrete,each=4) + rep(c(-offset,-offset,offset,offset),length(time.discrete))
km.points=c(0,km.points)

km.fun = list()

for(i in 1:3){
  surv.tmp = survival.fun[[i]]
  km.fun[[i]]=c(rep(surv.tmp[1:(length(surv.tmp)-1)],each=2),surv.tmp[length(surv.tmp)])
}

index.g = 3+2*(1:length(time.discrete)-1)
index.d = 2+2*(1:length(time.discrete)-1)

# Survival function plots for different combinations of germination and mortality
plot(NA,NA,
     ylim=c(0,1),xlim=c(0,max(time.discrete)), pch=19,axes=FALSE,
     ylab= "P(seed remains in seed bank)",xlab="Time (t)")
rect(xleft=time.discrete-.05,ybottom=-.25,xright=time.discrete+.15,ytop=1.25,border=0,col='gray95')

# Lines similar to Kaplan-Meier plot
for(i in 1:3){
  lines(km.points,km.fun[[i]],lty=plot.lines[i],col=plot.colors[i])
  points(time.points[index.g],survival.fun[[i]][index.g],pch=21,bg='white',col=plot.colors[i])
  points(time.points[index.d],survival.fun[[i]][index.d],pch=19,col=plot.colors[i])
}

mtext("Time (t)", side = 1,line=2.2)
axis(1L);axis(2L);

par(mfrow=c(2,1))
par(mar=c(1,4,0,2),
    oma=c(4,0,2,0))

plot(NA,NA, 
     ylab = "Unconditional emergence",
     xlab = "", xaxt='n',axes=FALSE,
     xlim=c(0,max(time.discrete)),ylim=c(0,.5),
     pch=19)

for(i in 1:3){
  points(time.discrete,pmf.g[[i]],
         col=plot.colors[i],lty=plot.lines[i], pch = 19,
         type='b')
}

axis(1L,labels=FALSE);axis(2L);

plot(NA,NA, 
     xlab = "Time (t)", ylab = "Unconditional mortality",
     xlim=c(0,max(time.discrete)), ylim=c(0,.5),
     axes=FALSE,
     pch=19)

for(i in 1:3){
  points(time.discrete,pmf.m[[i]],
         col=plot.colors[i],lty=plot.lines[i], pch = 19,
         type='b')
}

axis(1L);axis(2L);

dev.off()