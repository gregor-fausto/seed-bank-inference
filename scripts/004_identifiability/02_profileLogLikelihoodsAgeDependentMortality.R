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
primaryDirectory <- "../"
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
  p.m1 = parms[2];
  p.m2 = parms[3];
  n = dat$n.s_obs;
  y = dat$y.g_obs;
  t = dat$t_obs
  
  index=(1:length(t))[t==1]
  prob = p.g*(1-p.m1)
  tmp = -sum(dbinom(y[index], n[index], prob, log=TRUE))
  
  index=(1:length(t))[t==2]
  prob = p.g*(1-p.m2)*(1-p.g)*(1-p.m1)
  tmp2 = -sum(dbinom(y[index], n[index], prob, log=TRUE))
  
  tmp = tmp+tmp2
  return(tmp)
}


# likSeedAddCC3 = function(parms, dat){
#   p.g = parms[1];
#   p.m = parms[2];
#   n = dat$n.s_obs;
#   y = dat$y.g_obs;
#   t = dat$t_obs
#   
#   index=(1:length(t))[t==1]
#   prob = p.g*(1-p.m)
#   tmp = -sum(dbinom(y[index], n[index], prob, log=TRUE))
#   
#   index=(1:length(t))[t==2]
#   prob = p.g*(1-p.m)*(1-p.g)*(1-p.m)
#   tmp2 = -sum(dbinom(y[index], n[index], prob, log=TRUE))
#   
#   index=(1:length(t))[t==3]
#   prob = p.g*(1-p.m)*(1-p.g)*(1-p.m)*(1-p.g)*(1-p.m)
#   tmp3 = -sum(dbinom(y[index], n[index], prob, log=TRUE))
#   
#   tmp = tmp+tmp2+tmp3
#   return(tmp)
# }


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
  p.m1 = parms[2];
  p.m2 = parms[3];
  n = dat$n.s_obs;
  y = dat$y.s_obs;
  n.g = dat$n.g_obs;
  y.g = dat$y.g_obs
  t = dat$t_obs
  
  index=(1:length(t))[t==1]
  prob = (1-p.m1)
  tmp = -sum(dbinom(y.g[index], n.g[index], p.g, log=TRUE)+dbinom(y[index], n[index], prob))
  
  index=(1:length(t))[t==2]
  prob = (1-p.m2)*(1-p.g)*(1-p.m1)
  tmp2 = -sum(dbinom(y.g[index], n.g[index], p.g, log=TRUE)+dbinom(y[index], n[index], prob))
  
  tmp = tmp+tmp2
  return(tmp)
}

# 
# likSeedBagCC3 = function(parms, dat){
#   p.g = parms[1];
#   p.m = parms[2];
#   n = dat$n.s_obs;
#   y = dat$y.s_obs;
#   n.g = dat$n.g_obs;
#   y.g = dat$y.g_obs
#   t = dat$t_obs
#   
#   index=(1:length(t))[t==1]
#   prob = (1-p.m)
#   tmp = -sum(dbinom(y.g[index], n.g[index], p.g, log=TRUE)+dbinom(y[index], n[index], prob))
#   
#   index=(1:length(t))[t==2]
#   prob = (1-p.m)*(1-p.g)*(1-p.m)
#   tmp2 = -sum(dbinom(y.g[index], n.g[index], p.g, log=TRUE)+dbinom(y[index], n[index], prob))
#   
#   index=(1:length(t))[t==3]
#   prob = (1-p.m)*(1-p.g)*(1-p.m)*(1-p.g)*(1-p.m)
#   tmp3 = -sum(dbinom(y.g[index], n.g[index], p.g, log=TRUE)+dbinom(y[index], n[index], prob))
#   
#   tmp = tmp+tmp2+tmp3
#   return(tmp)
# }

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

# set font sizes
pt12 = 1
pt10 = 10/12
pt9 = 9/12
pt8 = 8/12
pt7 = 7/12
pt6 = 6/12


# - ++1 year ----

# x = y = seq(0.08,.12,by=0.001)
# 
# resmat = matrix(nrow=length(x),ncol=length(y))
# 
# for(i in 1:length(x)){
#   for(j in 1:length(y)){
#     resmat[i,j] = likSeedBagCC1(parms=c(x[i],y[j]),dat=data)
#   }
# }
# 
# min(resmat)
# 
# prof_1 = apply(-resmat,1,max)
# prof_2 = apply(-resmat,2,max)
# 
# plot(x,prof_1,type='l')
# plot(y,prof_2,type='l')
# 
# 
# matrix.image(t(resmat),x=x,y=y,do.contour=TRUE,bw=TRUE, nl = 8,extraLevel = c(360,370))
# # contour(x, y, resmat,nlevels=15)
# # contour(x, y, resmat,levels=c(seq(360,370,by=2),seq(370,390,by=10)),add=TRUE)
# points(.1,.1,pch=4,cex=2)
# 
# mtext(expression(p[m]),side=2,line=1.25,cex=pt9);
# mtext("Seed bag burial",cex=pt10,adj=.5,outer=FALSE,side=2,line=2.25)
# 
# mtext(paste("# observation times = ",1),cex=pt10,adj=.5)

# - ++2 years ----

panelLabs =  c("G.","H.","I.","J.","K.","L.")

tiff(filename=paste0("~/Dropbox/chapter-4/analysis/products/figures/identifiability-joint-variable-likelihood.tif"),
     height=3,width=6,units="in",res=300,compression="lzw",pointsize=12)

par(mfrow=c(2,3),mar=c(1.5,1.5,1,.5),oma=c(1,2,0,0)+.1,mgp=c(3,.45,0))


x = y = z = seq(0,.2,by=0.001)

resmat = array(dim=c(length(x),length(y),length(z)))

for(i in 1:length(x)){
  for(j in 1:length(y)){
    for(k in 1:length(z)){
      resmat[i,j,k] = likSeedBagCC2(parms=c(x[i],y[j],z[k]),dat=data)
    }
  }
}


# profile for pg
prof_1 = apply(resmat,1,min)

# profile for pm1
prof_2 = apply(resmat,2,min)

# profile for pm2
prof_3 = apply(resmat,3,min)

plot(x,prof_2,type='l',ylab="",xlab="",cex.axis=pt8);abline(v=0.1,lty=3)
mtext(side=3, line =0,paste0(panelLabs[1]), cex = pt8 , adj = 0)
mtext("Negative log-likelihood",side=2,line=1.25,cex=pt9);
mtext("Seed bag burial",cex=pt10,adj=.5,outer=FALSE,side=2,line=2.25)

plot(y,prof_3,type='l',ylab="",xlab="",cex.axis=pt8);abline(v=0.1,lty=3)
mtext(side=3, line =0,paste0(panelLabs[2]), cex = pt8 , adj = 0)

plot(y,prof_1,type='l',ylab="",xlab="",cex.axis=pt8);abline(v=0.1,lty=3)
mtext(side=3, line =0,paste0(panelLabs[3]), cex = pt8 , adj = 0)

# 
# min(resmat)
# matrix.image(t(resmat[,,40]),x=x,y=y,do.contour=TRUE,bw=TRUE, nl = 6,extraLevel = c(710,720))
# #contour(x, y, resmat,nlevels=15)
# #contour(x, y, resmat,levels=c(seq(700,714,by=2),seq(720,780,by=20)),add=TRUE)
# points(.1,.1,pch=4,cex=2)


# 2 years
x = y = z = seq(0,1,by=0.01)

resmat = array(dim=c(length(x),length(y),length(z)))

for(i in 1:length(x)){
  for(j in 1:length(y)){
    for(k in 1:length(z)){
      resmat[i,j,k] = likSeedAddCC2(parms=c(x[i],y[j],z[k]),dat=data)
    }
  }
}


# profile for pg
prof_1 = apply(resmat,1,min)

# profile for pm1
prof_2 = apply(resmat,2,min)

# profile for pm2
prof_3 = apply(resmat,3,min)

plot(x,prof_2,type='l',ylab="",xlab="",cex.axis=pt8);abline(v=0.1,lty=3)
mtext(side=3, line =0,paste0(panelLabs[4]), cex = pt8 , adj = 0)
mtext("Negative log-likelihood",side=2,line=1.25,cex=pt9);
mtext("Seed addition",cex=pt10,adj=.5,outer=FALSE,side=2,line=2.25)
mtext(expression(p[m1]),side=1,line=1.4,cex=pt9);

plot(y,prof_3,type='l',ylab="",xlab="",cex.axis=pt8);abline(v=0.1,lty=3)
mtext(side=3, line =0,paste0(panelLabs[5]), cex = pt8 , adj = 0)
mtext(expression(p[m2]),side=1,line=1.4,cex=pt9);

plot(y,prof_1,type='l',ylab="",xlab="",cex.axis=pt8);abline(v=0.1,lty=3)
mtext(side=3, line =0,paste0(panelLabs[6]), cex = pt8 , adj = 0)
mtext(expression(p[g]),side=1,line=1.4,cex=pt9);

dev.off()

# - +3 year ----
# 
# resmat = matrix(nrow=length(x),ncol=length(y))
# 
# for(i in 1:length(x)){
#   for(j in 1:length(y)){
#     resmat[i,j] = likSeedBagCC3(parms=c(x[i],y[j]),dat=data)
#   }
# }
# 
# min(resmat)
# matrix.image(t(resmat),x=x,y=y,do.contour=TRUE,bw=TRUE, nl = 6,extraLevel = c(1045))
# #contour(x, y, resmat,nlevels=15)
# #contour(x, y, resmat,levels=c(seq(1045,1060,by=5),seq(1070,1090,by=10)),add=TRUE)
# points(.1,.1,pch=4,cex=2)
# 
# mtext(paste("# observation times = ",3),cex=pt10,adj=.5)


# - +Seed addition experiment ----
subset <- dataFull$t_obs<=3
data <- cbind(dataFull$y.s_obs[subset],
              dataFull$y.g_obs[subset],
              dataFull$n.g_obs[subset],
              dataFull$n.s_obs[subset],
              dataFull$t_obs[subset])
colnames(data) = c("y.s_obs","y.g_obs","n.g_obs","n.s_obs","t_obs")
data = data.frame(data)


x = y = seq(0,.1,by=0.01)

resmat = matrix(nrow=length(x),ncol=length(y))

for(i in 1:length(x)){
  for(j in 1:length(y)){
    resmat[i,j] = likSeedAddCC1(parms=c(x[i],y[j]),dat=data)
  }
}

prof_1 = apply(-resmat,1,max)
prof_2 = apply(-resmat,2,max)

plot(x,prof_1,type='l')
plot(y,prof_2,type='l')

min(resmat)
matrix.image(t(resmat),x=x,y=y,do.contour=TRUE,bw=TRUE, nl = 6,extraLevel = c(400,2000,20000))
#contour(x, y, resmat,levels=c(380,500,1000,5000,20000))
points(.1,.1,pch=4,cex=2)

mtext(expression(p[g]),side=1,line=1.5,cex=pt9)
mtext(expression(p[m]),side=2,line=1.5,cex=pt9);
mtext("Seed addition",cex=pt10,adj=.5,outer=FALSE,side=2,line=2.5)

# 2 years
x = y = z = seq(0,1,by=0.01)

resmat = array(dim=c(length(x),length(y),length(z)))

for(i in 1:length(x)){
  for(j in 1:length(y)){
    for(k in 1:length(z)){
      resmat[i,j,k] = likSeedAddCC2(parms=c(x[i],y[j],z[k]),dat=data)
    }
  }
}


prof_1 = apply(-resmat,1,max)
prof_2 = apply(-resmat,2,max)
prof_3 = apply(-resmat,3,max)

plot(x,prof_1,type='l')
plot(y,prof_2,type='l')
plot(y,prof_3,type='l')

min(resmat)
matrix.image(t(resmat),x=x,y=y,do.contour=TRUE,bw=TRUE, nl = 6,extraLevel = c(720,800,2000,20000))
#contour(x, y, resmat)
#contour(x,y,resmat,add=TRUE,levels=c(750,900,2000,3000,5000),col='red')
points(.1,.1,pch=4,cex=2)

mtext(expression(p[g]),side=1,line=1.5,cex=pt9)


x = y = seq(.06,.14,by=0.001)

resmat = matrix(nrow=length(x),ncol=length(y))

for(i in 1:length(x)){
  for(j in 1:length(y)){
    resmat[i,j] = likSeedAddCC3(parms=c(x[i],y[j]),dat=data)
  }
}

min(resmat)
matrix.image(t(resmat),x=x,y=y,do.contour=TRUE,bw=TRUE, nl = 3,extraLevel = c(1060,1100,2000,20000,50000))
#contour(x, y, resmat,xlim=c(0,.2),ylim=c(0,.2))
#contour(x,y,resmat,add=TRUE,levels=c(1200,1500,2000,3000,5000),col='red')
points(.1,.1,pch=4,cex=2)
mtext(expression(p[g]),side=1,line=1.5,cex=pt9)

dev.off()
# - +2 year ----

# - +Seed bag burial experiment ----

