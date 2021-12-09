####
####
# Study parameter identifiability in non-parametric models
# this is for a constant mortality constant germination model
# with different number of time points
####
####

# - Libraries ----
library(MCMCvis)


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

# - +Extract posterior----
f.extract = function(file=x,params=c("lambda")){
  samples.rjags = readRDS(file);
  mat <- MCMCchains(samples.rjags,params=params);
  return(mat)
}


# - Extract posteriors ----

# - +Seed bag burial experiment ----

params.Np.bags.all = list()
for(i in 1:3){
  # get filenames
  fileNamesBags=list.files(paste0("outputs/002_statisticalModelFitting/01_identifiability/posteriors-NpVmCg-",i,"/seedBagBurial/"))
  
  # construct parameter table from file names
  parm.table.full.bags=data.frame(do.call(rbind,lapply(fileNamesBags,f)))
  
  # subset parameter table
  index.bags=parm.table.full.bags$p.m==0.1 & parm.table.full.bags$p.g==0.1 & parm.table.full.bags$n.bags == 100
  parm.table.bags = parm.table.full.bags[index.bags,]
  
  # list of full file paths
  fits.Np.bags = paste0("outputs/002_statisticalModelFitting/01_identifiability/posteriors-NpVmCg-",i,"/seedBagBurial/",fileNamesBags[index.bags])
  
  # extract parameters for mortality and germination
  params.Np.bags = lapply(fits.Np.bags,f.extract,params=c("p.m","p.g"))
  
  params.Np.bags.all[[i]] <- params.Np.bags
}

# - +Seed addition experiment ----

params.Np.add.all = list()
for(i in 1:3){
  # get filenames
  fileNamesAddition=list.files(paste0("outputs/002_statisticalModelFitting/01_identifiability/posteriors-NpVmCg-",i,"/seedAddition/"))
  
  # construct parameter table from file names
  parm.table.full.add=data.frame(do.call(rbind,lapply(fileNamesAddition,f)))
  
  # subset parameter table
  index.add=parm.table.full.add$p.m==0.1 & parm.table.full.add$p.g==0.1 & parm.table.full.bags$n.bags == 100
  parm.table.add = parm.table.full.add[index.add,]
  
  # list of full file paths
  fits.Np.add = paste0("outputs/002_statisticalModelFitting/01_identifiability/posteriors-NpVmCg-",i,"/seedAddition/",fileNamesAddition[index.add])
  
  # extract parameters for mortality and germination
  params.Np.add = lapply(fits.Np.add,f.extract,params=c("p.m","p.g"))
  
  params.Np.add.all[[i]] <- params.Np.add
}

# - Calculate and summarize posterior error ----

# - Plots ----
# set font sizes
pt12 = 1
pt10 = 10/12
pt9 = 9/12
pt8 = 8/12
pt7 = 7/12
pt6 = 6/12

# - +Prior-posterior overlap ----

# - ++Mortality ----

tiff(filename=paste0("products/figures/prior-posterior-overlap-variable-mortality.tif"),
     height=3,width=6,units="in",res=300,compression="lzw",pointsize=12)

f.pp = function(index=1){
  n = length(fits.Np.bags)
  pp.add = pp.bags = matrix(NA,nrow=n,ncol=3)
  for(i in 1:3){
    for(j in 1:n){
      prior = rbeta(1500,1,1)
      
      # seed bags: mortality
      post = params.Np.bags.all[[i]][[j]][,i]
      pp <- list(prior, post)
      pp.bags[j,i]<-overlapping::overlap(pp)$OV[[1]]
      
      # seed addition: mortality
      post = params.Np.add.all[[i]][[j]][,i]
      pp <- list(prior, post)
      pp.add[j,i]<-overlapping::overlap(pp)$OV[[1]]
    }
  }
  return(list(pp.bags,pp.add))
}

par(mfrow = c(1,2),
    oma = c(1.5,2.5,0,0) + 0.1,
    mar = c(1,0,1,0.25) + 0.1)

plot(NA,
     xlim=c(.5,3.5),ylim=c(0,.75),
     col=c('black'),pch=c(16),
     xlab="",
     ylab="",
     xaxt='n',yaxt='n')


for(timepoint in 1:3){
  
  x = f.pp(index=timepoint)
  pp.bags = x[[1]]; pp.add = x[[2]]
  color.palette <- colorRampPalette(c("orange","purple"))(5)
  
  abline(h=0.35,lty='dotted')
  
  if(timepoint==1){offset=c(0,-.125,-.25)} else if(timepoint==2){offset=c(0,.125,0)} else { offset=c(0,0,.25)}
  pchVec = c(16,15,17)
  pchVec2 = c(21, 22, 24)
  for(i in (1:3)[(1:3)>=timepoint]){
    points(x=rep(i,1)+offset[i],y=pp.bags[,i],     
           col="black",pch=pchVec[timepoint])
    points(x=rep(i,1)+offset[i],y=pp.add[,i],
           col="black",pch=pchVec2[timepoint])
  }
  
  
}
axis(1, at = c(0,1,2,3), cex.axis = pt10,padj=-.75)
axis(2, at = seq(0,.7,by=.1), las = 1, cex.axis = pt10,hadj=.75)

L = legend("topright",
           title="Parameter",ncol = 2,
           legend=rep(NA,6),horiz=FALSE, bty='n', 
           pch=c(16,15,17, 21, 22, 24),col=c("black", "black"), inset=c(.05,0),
           cex=pt8,bg='gray',box.lty=0, x.intersp = .25)

legend(x = L$rect$left*1.125, y = L$rect$top*.93, legend = c(expression(p[m1]),expression(p[m2]),expression(p[m3])), 
       ncol=1, bg = NA, bty = 'n',cex=pt8,inset=c(.05,0))


box()
mtext( expression(Prior-posterior ~ overlap ), side=2, line = 1.8, cex =pt10,adj=.5)
mtext( "Years of observations", side=1, line = .3, cex =pt10, outer = TRUE)

mtext( "C.", side=3, line = 0, cex =pt10,adj=0)


# - ++Germination ----
n = length(fits.Np.bags)
pp.add = pp.bags = matrix(NA,nrow=n,ncol=5)
for(i in 1:3){
  for(j in 1:n){
    prior = rbeta(1500,1,1)
    
    # seed bags: mortality
    post = params.Np.bags.all[[i]][[j]][,6]
    pp <- list(prior, post)
    pp.bags[j,i]<-overlapping::overlap(pp)$OV[[1]]
    
    # seed addition: mortality
    post = params.Np.add.all[[i]][[j]][,6]
    pp <- list(prior, post)
    pp.add[j,i]<-overlapping::overlap(pp)$OV[[1]]
  }
}



plot(NA,
     xlim=c(.5,3.5),ylim=c(0,.75),
     col=c('black'),pch=c(16),
     xlab="",
     ylab="",
     xaxt='n',yaxt='n')

for(timepoint in 1:3){
  
  
  
  x = f.pp(index=timepoint)
  pp.bags = x[[1]]; pp.add = x[[2]]
  color.palette <- colorRampPalette(c("orange","purple"))(5)
  
  abline(h=0.35,lty='dotted')
  
  if(timepoint==1){offset=c(0,-.125,-.25)} else if(timepoint==2){offset=c(0,.125,0)} else { offset=c(0,0,.25)}
  pchVec = c(16,15,17)
  pchVec2 = c(21, 22, 24)
  for(i in (1:3)[(1:3)>=timepoint]){
    points(x=rep(i,1)+offset[i],y=pp.bags[,i],     
           col="black",pch=pchVec[timepoint])
    points(x=rep(i,1)+offset[i],y=pp.add[,i],
           col="black",pch=pchVec2[timepoint])
  }
  
}
axis(1, at = c(0,1,2,3), cex.axis = pt10,padj=-.75)

L = legend("topright",
           title="Parameter",ncol = 2,
           legend=rep(NA,6),horiz=FALSE, bty='n', 
           pch=c(16,15,17, 21, 22, 24),col=c("black", "black"), inset=c(.05,0),
           cex=pt8,bg='gray',box.lty=0, x.intersp = .25)

legend(x = L$rect$left*1.125, y = L$rect$top*.93, legend = c(expression(p[g1]),expression(p[g2]),expression(p[g3])), 
       ncol=1, bg = NA, bty = 'n',cex=pt8,inset=c(.05,0))


box()
mtext( "D.", side=3, line = 0, cex =pt10,adj=0)


dev.off()

# - +Correlations in joint posterior ----


tiff(filename=paste0("products/figures/identifiability-joint-variable-1.tif"),
     height=3,width=2,units="in",res=300,compression="lzw",pointsize=12)

par(mfrow=c(2,1),mar=c(.25,.25,.6,0),oma=c(1,2.5,.5,.25)+.1,mgp=c(3,.45,0))

i=1
# seed bags: mortality
post.bags = params.Np.bags.all[[i]][[1]][,c(1:i,6)]

plot(post.bags[,2],post.bags[,1],pch=16,ylab="",xlab="",cex=.5,
     xlim=c(.097,.105),ylim=c(.097,.105),tick=TRUE,
     axes=FALSE)

axis(2,cex.axis=pt6)

mtext(expression(p[m1]),side=2,line=1.25,cex=pt9);
mtext("Seed bag burial",cex=pt10,adj=.5,outer=FALSE,side=2,line=2)

mtext(adj=.7,paste("1 year of observations"),cex=pt8,outer=TRUE,line=-.5)
mtext(adj=0,"G.",cex=pt10,outer=FALSE)

box()

# seed addition: mortality
post.add = params.Np.add.all[[i]][[1]][,c(1:i,6)]

plot(post.add[,2],post.add[,1],pch=21,cex=.5,xlab="",ylab="",
     xlim=c(0.1,1),ylim=c(.1,1),tick=TRUE,
     axes=FALSE)

axis(2,cex.axis=pt6)

mtext(expression(p[m1]),side=2,line=1.25,cex=pt9)
mtext(expression(p[g]),side=1,line=0,cex=pt9);
mtext("Seed addition",cex=pt10,adj=.5,outer=FALSE,side=2,line=2)
mtext(adj=0,"J.",cex=pt10,outer=FALSE)

box()

dev.off()

tiff(filename=paste0("products/figures/identifiability-joint-variable-2.tif"),
     height=3,width=2,units="in",res=300,compression="lzw",pointsize=12)

# age 2

par(mfrow=c(4,2),mar=c(.25,.25,0,0),oma=c(1,1.75,.8,0)+.1,mgp=c(3,.45,0))


i=2
# seed bags: mortality
post.bags = params.Np.bags.all[[i]][[1]][,c(1:i,6)]

plot(post.bags[,2],post.bags[,1],pch=16,cex=.4,ylab="",xlab="",
     xlim=c(.092,.106),ylim=c(.092,.106),
     xaxt='n',yaxt='n')

axis(2,cex.axis=pt6,tick=TRUE,outer=FALSE,at=seq(.09,.11,by=.005))
mtext(adj=0,"H.",cex=pt10,outer=FALSE)



mtext( expression(p[m1] ), side=2, line = 1, cex =pt10,adj=0.5,outer=FALSE)
mtext( expression(p[m2] ), side=1, line = 0.5, cex =pt10,adj=0.5,outer=FALSE)
box()

plot(post.bags[,3],post.bags[,1],pch=16,cex=.4,ylab="",xlab="",
     xlim=c(.092,.106),ylim=c(.092,.106),
     xaxt='n',yaxt='n')

plot(NA,xlim=c(0,1),ylim=c(0,1),axes=FALSE)

plot(post.bags[,3],post.bags[,2],pch=16,cex=.4,ylab="",xlab="",
     xlim=c(.092,.106),ylim=c(.092,.106),
     xaxt='n',yaxt='n')

axis(2,cex.axis=pt6,tick=TRUE,outer=FALSE,at=seq(.09,.11,by=.005))

mtext( expression(p[m2] ), side=2, line = 1.25, cex =pt10,adj=0.5,outer=FALSE)
mtext(outer=TRUE,"2 years of observations",side=3, cex = pt8)

# seed addition: mortality
post.add = params.Np.add.all[[i]][[1]][,c(1:i,6)]

plot(post.add[,2],post.add[,1],pch=21,cex=.4,ylab="",xlab="",
     xlim=c(0,.6),ylim=c(0,.6),
     xaxt='n',yaxt='n')
mtext(adj=0,"K.",cex=pt10,outer=FALSE)

axis(2,cex.axis=pt6,outer=FALSE,seq(0,.6,by=.2))

mtext( expression(p[m1] ), side=2, line = 1, cex =pt10,adj=0.5,outer=FALSE)
mtext( expression(p[m2] ), side=1, line = 0.5, cex =pt10,adj=0.5,outer=FALSE)

plot(post.add[,3],post.add[,1],pch=16,cex=.4,ylab="",xlab="",
     xlim=c(0,.6),ylim=c(0,.6),
     xaxt='n',yaxt='n')

plot(NA,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',axes=FALSE)

plot(post.add[,3],post.add[,2],pch=21,cex=.4,ylab="",xlab="",
     xlim=c(0,.2),ylim=c(0,.2),
     xaxt='n',yaxt='n')

axis(2,cex.axis=pt6,outer=FALSE,seq(0,.2,by=.1))

mtext( expression(p[m2] ), side=2, line = 1.25, cex =pt10,adj=0.5,outer=FALSE)
mtext( expression(p[g] ), side=1, line = 0.5, cex =pt10,adj=0.5,outer=FALSE)

dev.off()


tiff(filename=paste0("products/figures/identifiability-joint-variable-3.tif"),
     height=3,width=2,units="in",res=300,compression="lzw",pointsize=12)
# age 3

par(mfrow=c(6,3),mar=c(.25,.25,.25,0),oma=c(1,1.75,.8,0)+.1,mgp=c(3,.45,0))

i=3
# seed bags: mortality
post.bags = params.Np.bags.all[[i]][[1]][,c(1:i,6)]

row1= c(.09,.115)

plot(post.bags[,2],post.bags[,1],pch=16,cex=.4,ylab="",xlab="",
     xlim=row1,ylim=row1,
     xaxt='n',yaxt='n')

axis(2,cex.axis=pt6,tick=TRUE,outer=FALSE,at=seq(.09,.12,by=.01))

mtext( expression(p[m1] ), side=2, line = 1, cex =pt10,adj=0.5,outer=FALSE)
mtext( expression(p[m2] ), side=1, line = .1, cex =pt10,adj=0.5,outer=FALSE)
mtext(adj=0,"I.",cex=pt10,outer=FALSE)


plot(post.bags[,3],post.bags[,1],pch=16,cex=.4,ylab="",xlab="",
     xlim=row1,ylim=row1,
     xaxt='n',yaxt='n')

plot(post.bags[,4],post.bags[,1],pch=16,cex=.4,ylab="",xlab="",
     xlim=row1,ylim=row1,
     xaxt='n',yaxt='n')

plot(NA,xlim=c(0,1),ylim=c(0,1),axes=FALSE)

mtext(outer=TRUE,"3 years of observations",side=3, cex = pt8)

row2 = c(.09,.115)


plot(post.bags[,3],post.bags[,1],pch=16,cex=.4,ylab="",xlab="",
     xlim=row2,ylim=row2,
     xaxt='n',yaxt='n')

axis(2,cex.axis=pt6,tick=TRUE,outer=FALSE,at=seq(.09,.12,by=.01))

mtext( expression(p[m2] ), side=2, line = 1, cex =pt10,adj=0.5,outer=FALSE)
mtext( expression(p[m3] ), side=1, line = .1, cex =pt10,adj=0.5,outer=FALSE)

plot(post.bags[,4],post.bags[,2],pch=16,cex=.4,ylab="",xlab="",
     xlim=row2,ylim=row2,
     xaxt='n',yaxt='n')

plot(NA,xlim=c(0,1),ylim=c(0,1),axes=FALSE)
plot(NA,xlim=c(0,1),ylim=c(0,1),axes=FALSE)

row3=c(.09,.115)

plot(post.bags[,4],post.bags[,3],pch=16,cex=.4,ylab="",xlab="",
     xlim=row3,ylim=row3,
     xaxt='n',yaxt='n')

axis(2,cex.axis=pt6,tick=TRUE,outer=FALSE,at=seq(.09,.12,by=.01))


mtext( expression(p[m3] ), side=2, line = 1, cex =pt10,adj=0.5,outer=FALSE)

# seed addition: mortality
post.add = params.Np.add.all[[i]][[1]][,c(1:i,6)]

row1=c(0,.55)

plot(post.add[,2],post.add[,1],pch=21,cex=.4,ylab="",xlab="",
     xlim=row1,ylim=row1,
     xaxt='n',yaxt='n')

axis(2,cex.axis=pt6,tick=TRUE,outer=FALSE,at=c(0,.2,.4))

mtext( expression(p[m1] ), side=2, line = 1, cex =pt10,adj=0.5,outer=FALSE)
mtext( expression(p[m2] ), side=1, line = .1, cex =pt10,adj=0.5,outer=FALSE)
mtext(adj=0,"L.",cex=pt10,outer=FALSE)

plot(post.add[,3],post.add[,1],pch=21,cex=.4,ylab="",xlab="",
     xlim=row1,ylim=row1,
     xaxt='n',yaxt='n')

plot(post.add[,4],post.add[,1],pch=21,cex=.4,ylab="",xlab="",
     xlim=row1,ylim=row1,
     xaxt='n',yaxt='n')

plot(NA,xlim=c(0,1),ylim=c(0,1),axes=FALSE)

row2=c(0,.2)

plot(post.add[,3],post.add[,2],
     pch=21,cex=.4,ylab="",xlab="",
     xlim=row2,ylim=row2,
     xaxt='n',yaxt='n')

axis(2,cex.axis=pt6,tick=TRUE,outer=FALSE,at=seq(0,.2,by=.1))

mtext( expression(p[m2] ), side=2, line = 1, cex =pt10,adj=0.5,outer=FALSE)
mtext( expression(p[m3] ), side=1, line = .1, cex =pt10,adj=0.5,outer=FALSE)

plot(post.add[,4],post.add[,2],pch=21,cex=.4,ylab="",xlab="",
     xlim=row2,ylim=row2,
     xaxt='n',yaxt='n')

plot(NA,xlim=c(0,1),ylim=c(0,1),axes=FALSE)
plot(NA,xlim=c(0,1),ylim=c(0,1),axes=FALSE)

row3 = c(0,.2)
plot(post.add[,4],post.add[,3],
     pch=21,cex=.4,ylab="",xlab="",
     xlim=row3,ylim=row3,
     xaxt='n',yaxt='n')

axis(2,cex.axis=pt6,tick=TRUE,outer=FALSE,at=seq(0,.2,by=.1))

mtext( expression(p[m3] ), side=2, line = 1, cex =pt10,adj=0.5,outer=FALSE)
mtext( expression(p[g] ), side=1, line = 0.5, cex =pt10,adj=0.5,outer=FALSE)

box()


dev.off()
