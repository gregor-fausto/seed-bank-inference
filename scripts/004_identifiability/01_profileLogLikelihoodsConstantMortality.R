####
####
# Study parameter identifiability in non-parametric models
# this is for a constant mortality constant germination model
# with different number of time points
####
####

# - Libraries ----
library(MCMCvis)

# - Directories ----
primaryDirectory <- ""
outDataDirectory <- paste0(primaryDirectory,"outputs/001_simulateObservations/01_identifiability/")
outPosteriorSamplesDirectory <- paste0(primaryDirectory,"outputs/002_statisticalModelFitting/01_identifiability/")

# - Functions ----

matrix.image=function(A, x=NULL, y=NULL, col=rainbow(100,start=0.67,end=0),
                      bw=FALSE, do.contour=FALSE, do.legend=FALSE,nl = 10, extraLevel = NULL ,...) {
  if(is.null(x)) x=1:ncol(A);
  if(is.null(y)) y=1:nrow(A); 
  nx=length(x); ny=length(y); 
  x1=c(1.5*x[1]-0.5*x[2],1.5*x[nx]-0.5*x[nx-1]); 
  y1=c(1.5*y[1]-0.5*y[2],1.5*y[ny]-0.5*y[ny-1]); 
  if(bw) col=grey( (200:0)/200 ); 
  # comment out this line to reverse the direction of the plot
  image(list(x=x,y=y,z=A),xlim=x1,ylim=(y1),col=col,cex.axis=pt8,cex.lab=pt8,bty="u",...);
  abline(v=range(x1)); abline(h=range(y1)); 
  if(do.contour) contour(x,y,A,nlevels=nl,labcex=pt6,add=TRUE);   
  contour(x,y,A,levels=extraLevel,labcex=pt6,add=TRUE);
}

# set font sizes
pt12 = 1
pt10 = 10/12
pt9 = 9/12
pt8 = 8/12
pt7 = 7/12
pt6 = 6/12


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
parametersIdentifiability <- parameterTable$p.m==0.1 & parameterTable$p.g==0.1 & parameterTable$n.bags %in% c(100);
index = (1:n)[parametersIdentifiability]

# - Fetch data ----

# - +Non-parametric, constant mortality
for(i in 1:length(index)){
  simulation.data  <-  readRDS(paste0(paste0(outDataDirectory),fileNames[index[i]]))
}

dataFull <- simulation.data[[2]]


# - Negative log-likelihood functions ----

# - +C/C seed addition ----
likSeedAddCC1 = function(parms, dat){
  p.g = parms[1];
  p.m = parms[2];
  n = dat$n.s_obs;
  y = dat$y.g_obs;
  t = dat$t_obs
  index=(1:length(t))[t==1]
  prob = p.g*(1-p.m)
  tmp = -sum(dbinom(y[index], n[index], prob, log=TRUE))
  return(tmp)
}

likSeedAddCC2 = function(parms, dat){
  p.g = parms[1];
  p.m = parms[2];
  n = dat$n.s_obs;
  y = dat$y.g_obs;
  t = dat$t_obs
  
  index=(1:length(t))[t==1]
  prob = p.g*(1-p.m)
  tmp = -sum(dbinom(y[index], n[index], prob, log=TRUE))
  
  index=(1:length(t))[t==2]
  prob = p.g*(1-p.m)*(1-p.g)*(1-p.m)
  tmp2 = -sum(dbinom(y[index], n[index], prob, log=TRUE))
  
  tmp = tmp+tmp2
  return(tmp)
}


likSeedAddCC3 = function(parms, dat){
  p.g = parms[1];
  p.m = parms[2];
  n = dat$n.s_obs;
  y = dat$y.g_obs;
  t = dat$t_obs
  
  index=(1:length(t))[t==1]
  prob = p.g*(1-p.m)
  tmp = -sum(dbinom(y[index], n[index], prob, log=TRUE))
  
  index=(1:length(t))[t==2]
  prob = p.g*(1-p.m)*(1-p.g)*(1-p.m)
  tmp2 = -sum(dbinom(y[index], n[index], prob, log=TRUE))
  
  index=(1:length(t))[t==3]
  prob = p.g*(1-p.m)*(1-p.g)*(1-p.m)*(1-p.g)*(1-p.m)
  tmp3 = -sum(dbinom(y[index], n[index], prob, log=TRUE))
  
  tmp = tmp+tmp2+tmp3
  return(tmp)
}


# - +C/C seed bag ----
likSeedBagCC1 = function(parms, dat){
  p.g = parms[1];
  p.m = parms[2];
  n = dat$n.s_obs;
  y = dat$y.s_obs;
  n.g = dat$n.g_obs;
  y.g = dat$y.g_obs
  t = dat$t_obs
  
  index=(1:length(t))[t==1]
  prob = (1-p.m)
  tmp = -sum(dbinom(y.g[index], n.g[index], p.g, log=TRUE)+dbinom(y[index], n[index], prob))
  return(tmp)
}

likSeedBagCC2 = function(parms, dat){
  p.g = parms[1];
  p.m = parms[2];
  n = dat$n.s_obs;
  y = dat$y.s_obs;
  n.g = dat$n.g_obs;
  y.g = dat$y.g_obs
  t = dat$t_obs
  
  index=(1:length(t))[t==1]
  prob = (1-p.m)
  tmp = -sum(dbinom(y.g[index], n.g[index], p.g, log=TRUE)+dbinom(y[index], n[index], prob))
  
  index=(1:length(t))[t==2]
  prob = (1-p.m)*(1-p.g)*(1-p.m)
  tmp2 = -sum(dbinom(y.g[index], n.g[index], p.g, log=TRUE)+dbinom(y[index], n[index], prob))
  
  tmp = tmp+tmp2
  return(tmp)
}


likSeedBagCC3 = function(parms, dat){
  p.g = parms[1];
  p.m = parms[2];
  n = dat$n.s_obs;
  y = dat$y.s_obs;
  n.g = dat$n.g_obs;
  y.g = dat$y.g_obs
  t = dat$t_obs
  
  index=(1:length(t))[t==1]
  prob = (1-p.m)
  tmp = -sum(dbinom(y.g[index], n.g[index], p.g, log=TRUE)+dbinom(y[index], n[index], prob))
  
  index=(1:length(t))[t==2]
  prob = (1-p.m)*(1-p.g)*(1-p.m)
  tmp2 = -sum(dbinom(y.g[index], n.g[index], p.g, log=TRUE)+dbinom(y[index], n[index], prob))
  
  index=(1:length(t))[t==3]
  prob = (1-p.m)*(1-p.g)*(1-p.m)*(1-p.g)*(1-p.m)
  tmp3 = -sum(dbinom(y.g[index], n.g[index], p.g, log=TRUE)+dbinom(y[index], n[index], prob))
  
  tmp = tmp+tmp2+tmp3
  return(tmp)
}

# - Simulate data ----

subset <- dataFull$t_obs<=3
data <- cbind(dataFull$y.s_obs[subset],
              dataFull$y.g_obs[subset],
              dataFull$n.g_obs[subset],
              dataFull$n.s_obs[subset],
              dataFull$t_obs[subset])
colnames(data) = c("y.s_obs","y.g_obs","n.g_obs","n.s_obs","t_obs")
data = data.frame(data)

# - +Seed bag burial experiment ----
panelLabs =  c("A.","B.","C.","D.","E.","F.")


tiff(filename=paste0("~/Dropbox/chapter-4/analysis/products/figures/identifiability-joint-constant-likelihood.tif"),
     height=3,width=6,units="in",res=300,compression="lzw",pointsize=12)

par(mfrow=c(2,3),mar=c(1.5,1.5,1,.5),oma=c(1,2,1,0)+.1,mgp=c(3,.45,0))

# - ++1 year ----

x = y = seq(0.08,.12,by=0.0001)

resmat = matrix(nrow=length(x),ncol=length(y))

for(i in 1:length(x)){
  for(j in 1:length(y)){
    resmat[i,j] = likSeedBagCC1(parms=c(x[i],y[j]),dat=data)
  }
}

matrix.image(resmat,x=x,y=y,do.contour=TRUE,bw=TRUE, nl = 8,extraLevel = c(360,370))
points(.1,.1,pch=4,cex=2)

mtext(expression(p[m]),side=2,line=1.25,cex=pt9);
mtext("Seed bag burial",cex=pt10,adj=.5,outer=FALSE,side=2,line=2.25)

mtext(side=3, line =0,paste0(panelLabs[1]), cex = pt8 , adj = 0)
mtext(side=3, line =0,paste0("1 year of observations"), cex = pt8 , adj = .5)

# - ++2 year ----

resmat = matrix(nrow=length(x),ncol=length(y))

for(i in 1:length(x)){
  for(j in 1:length(y)){
    resmat[i,j] = likSeedBagCC2(parms=c(x[i],y[j]),dat=data)
  }
}

matrix.image(resmat,x=x,y=y,do.contour=TRUE,bw=TRUE, nl = 6,extraLevel = c(710,720))
points(.1,.1,pch=4,cex=2)

mtext(side=3, line =0,paste0(panelLabs[2]), cex = pt8 , adj = 0)
mtext(side=3, line =0,paste0("2 years of observations"), cex = pt8 , adj = .5)

# - ++3 year ----

resmat = matrix(nrow=length(x),ncol=length(y))

for(i in 1:length(x)){
  for(j in 1:length(y)){
    resmat[i,j] = likSeedBagCC3(parms=c(x[i],y[j]),dat=data)
  }
}

matrix.image(resmat,x=x,y=y,do.contour=TRUE,bw=TRUE, nl = 6,extraLevel = c(1045))
points(.1,.1,pch=4,cex=2)

mtext(side=3, line =0,paste0(panelLabs[3]), cex = pt8 , adj = 0)
mtext(side=3, line =0,paste0("3 years of observations"), cex = pt8 , adj = .5)


# - +Seed addition experiment ----
subset <- dataFull$t_obs<=3
data <- cbind(dataFull$y.s_obs[subset],
              dataFull$y.g_obs[subset],
              dataFull$n.g_obs[subset],
              dataFull$n.s_obs[subset],
              dataFull$t_obs[subset])
colnames(data) = c("y.s_obs","y.g_obs","n.g_obs","n.s_obs","t_obs")
data = data.frame(data)

# - ++1 year ----

x = y = seq(0,1,by=0.01)

resmat = matrix(nrow=length(x),ncol=length(y))

for(i in 1:length(x)){
  for(j in 1:length(y)){
    resmat[i,j] = likSeedAddCC1(parms=c(x[i],y[j]),dat=data)
  }
}

matrix.image(resmat,x=x,y=y,do.contour=TRUE,bw=TRUE, nl = 6,extraLevel = c(400,2000,20000))
points(.1,.1,pch=4,cex=2)

mtext(expression(p[g]),side=1,line=1.5,cex=pt9)
mtext(expression(p[m]),side=2,line=1.5,cex=pt9);

mtext("Seed addition",cex=pt10,adj=.5,outer=FALSE,side=2,line=2.5)
mtext(side=3, line =0,paste0(panelLabs[4]), cex = pt8 , adj = 0)

# - ++2 year ----
x = y = seq(0.06,.14,by=0.0001)

resmat = matrix(nrow=length(x),ncol=length(y))

for(i in 1:length(x)){
  for(j in 1:length(y)){
    resmat[i,j] = likSeedAddCC2(parms=c(x[i],y[j]),dat=data)
  }
}

matrix.image(resmat,x=x,y=y,do.contour=TRUE,bw=TRUE, nl = 6,extraLevel = c(720,800,2000,20000))
points(.1,.1,pch=4,cex=2)

mtext(expression(p[g]),side=1,line=1.5,cex=pt9)
mtext(side=3, line =0,paste0(panelLabs[5]), cex = pt8 , adj = 0)

# - ++3 year ----

x = y = seq(.06,.14,by=0.0001)

resmat = matrix(nrow=length(x),ncol=length(y))

for(i in 1:length(x)){
  for(j in 1:length(y)){
    resmat[i,j] = likSeedAddCC3(parms=c(x[i],y[j]),dat=data)
  }
}

matrix.image(resmat,x=x,y=y,do.contour=TRUE,bw=TRUE, nl = 3,extraLevel = c(1060,1100,2000,20000,50000))
points(.1,.1,pch=4,cex=2)

mtext(expression(p[g]),side=1,line=1.5,cex=pt9)
mtext(side=3, line =0,paste0(panelLabs[6]), cex = pt8 , adj = 0)

dev.off()
