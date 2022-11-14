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

primaryDirectory <- ""

# - +Seed bag burial experiment ----

params.Np.bags.all = list()
for(i in 1:3){
  # get filenames
  fileNamesBags=list.files(paste0(primaryDirectory,"outputs/002_statisticalModelFitting/01_identifiability/posteriors-NpCmCg-",i,"/seedBagBurial/"))

  # construct parameter table from file names
  parm.table.full.bags=data.frame(do.call(rbind,lapply(fileNamesBags,f)))

  # subset parameter table
  index.bags=parm.table.full.bags$p.m==0.1 & parm.table.full.bags$p.g==0.1
  parm.table.bags = parm.table.full.bags[index.bags,]

  # list of full file paths
  fits.Np.bags = paste0(primaryDirectory,"outputs/002_statisticalModelFitting/01_identifiability/posteriors-NpCmCg-",i,"/seedBagBurial/",fileNamesBags[index.bags])

  # extract parameters for mortality and germination
  params.Np.bags = lapply(fits.Np.bags,f.extract,params=c("p.m","p.g"))

  params.Np.bags.all[[i]] <- params.Np.bags
}

# - +Seed addition experiment ----

params.Np.add.all = list()
for(i in 1:3){
  # get filenames
  fileNamesAddition=list.files(paste0(primaryDirectory,"outputs/002_statisticalModelFitting/01_identifiability/posteriors-NpCmCg-",i,"/seedAddition/"))

  # construct parameter table from file names
  parm.table.full.add=data.frame(do.call(rbind,lapply(fileNamesAddition,f)))

  # subset parameter table
  index.add=parm.table.full.add$p.m==0.1 & parm.table.full.add$p.g==0.1
  parm.table.add = parm.table.full.add[index.add,]

  # list of full file paths
  fits.Np.add = paste0(primaryDirectory,"outputs/002_statisticalModelFitting/01_identifiability/posteriors-NpCmCg-",i,"/seedAddition/",fileNamesAddition[index.add])

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

# - +Prior-posterior overlap ----

# - ++Mortality ----

tiff(filename=paste0(primaryDirectory,"products/figures/prior-posterior-overlap-constant.tif"),
     height=3,width=6,units="in",res=300,compression="lzw",pointsize=12)
par(mfrow = c(1,2),
    oma = c(1.5,2.5,0,0) + 0.1,
    mar = c(1,0,1,0.25) + 0.1)

n = length(fits.Np.bags)
pp.add = pp.bags = matrix(NA,nrow=n,ncol=3)
for(i in 1:3){
  for(j in 1:n){
  prior = rbeta(1500,1,1)

  # seed bags: mortality
  post = params.Np.bags.all[[i]][[j]][,1]
  pp <- list(prior, post)
  pp.bags[j,i]<-overlapping::overlap(pp)$OV[[1]]

  # seed addition: mortality
  post = params.Np.add.all[[i]][[j]][,1]
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

axis(1, at = c(0,1,2,3), cex.axis = pt10,padj=-.75)
axis(2, at = seq(0,.7,by=.1), las = 1, cex.axis = pt10,hadj=.75)

mtext( "A.", side=3, line = 0, cex =pt10,adj=0)
mtext( expression(Prior-posterior ~ overlap ), side=2, line = 1.8, cex =pt10,adj=.5)

abline(h=0.35,lty='dotted')

for(i in 1:3){
  points(x=rep(i,1),y=pp.bags[,i],
         col="black",pch=c(16))
  points(x=rep(i,1),y=pp.add[,i],
         col="black",pch=c(21))
}



L = legend("topright",
           title="Parameter",ncol = 2,
           legend=rep(NA,2),horiz=FALSE, bty='n',
           pch=c(16, 21 ),col=c("black", "black"), inset=c(.05,0),
           cex=pt8,bg='gray',box.lty=0, x.intersp = .25)

legend(x = L$rect$left*1.125, y = L$rect$top*.93, legend = c(expression(p[m])),
       ncol=1, bg = NA, bty = 'n',cex=pt8,inset=c(.05,0))


box()


# - ++Germination ----
n = length(fits.Np.bags)
pp.add = pp.bags = matrix(NA,nrow=n,ncol=3)
for(i in 1:3){
  for(j in 1:n){
    prior = rbeta(1500,1,1)

    # seed bags: mortality
    post = params.Np.bags.all[[i]][[j]][,2]
    pp <- list(prior, post)
    pp.bags[j,i]<-overlapping::overlap(pp)$OV[[1]]

    # seed addition: mortality
    post = params.Np.add.all[[i]][[j]][,2]
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

axis(1, at = c(0,1,2,3), cex.axis = pt10,padj=-.75)

mtext( "B.", side=3, line = 0,cex=pt10,adj=0)

mtext( "Years of observations", side=1, line = .3, cex =pt10, outer = TRUE)

abline(h=0.35,lty='dotted')

for(i in 1:3){
  points(x=rep(i,1),y=pp.bags[,i],
         col="black",pch=c(16))
  points(x=rep(i,1),y=pp.add[,i],
         col="black",pch=c(21))
}


box()


L = legend("topright",
           title="Parameter",ncol = 2,
           legend=rep(NA,2),horiz=FALSE, bty='n',
           pch=c(16, 21 ),col=c("black", "black"), inset=c(.05,0),
           cex=pt8,bg='gray',box.lty=0, x.intersp = .25)

legend(x = L$rect$left*1.125, y = L$rect$top*.93, legend = c(expression(p[g])),
       ncol=1, bg = NA, bty = 'n',cex=pt8,inset=c(.05,0))


dev.off()


# # - +Correlations in joint posterior ----

tiff(filename=paste0(primaryDirectory,"products/figures/identifiability-joint-constant.tif"),
     height=3,width=6,units="in",res=300,compression="lzw",pointsize=12)

par(mfrow=c(2,3),mar=c(1.5,1.5,.7,0),oma=c(1,2,1,0)+.1,mgp=c(3,.45,0))

panelLabs=c("A.","B.","C.")

for(i in 1:3){
  # seed bags: mortality
  post.bags = params.Np.bags.all[[i]][[1]]

  plot(post.bags[,2],post.bags[,1],pch=16,cex=.5,
       ylab="",xlab="",xlim=c(0.08,.12),ylim=c(0.08,.12),
       cex.axis=pt8)

  if(i==1){mtext(expression(p[m]),side=2,line=1.25,cex=pt9);
    mtext("Seed bag burial",cex=pt10,adj=.5,outer=FALSE,side=2,line=2.5)}

  mtext(side=3, line =0,paste0(panelLabs[i]), cex = pt8 , adj = 0)

 if(i==1){
   mtext(paste(i," year of observations"),cex=pt8,adj=.5,padj=-.1)
 } else {
   mtext(paste(i," years of observations"),cex=pt8,adj=.5,padj=-.1)

   }

  box()

}

panelLabs=c("D.","E.","F.")

for(i in 1:3){
  # seed addition: mortality
  post.add = params.Np.add.all[[i]][[1]]

  if(i!=1){plot(post.add[,2],post.add[,1],pch=21,cex=.5,xlab="",ylab="",
                xlim=c(0,.25),ylim=c(0,.25),
                cex.axis=pt8)} else {
    plot(post.add[,2],post.add[,1],pch=21,cex=.5,xlab="",ylab="",cex.axis=pt8,
         xlim=c(0,1),ylim=c(0,1))
  }
   mtext(expression(p[g]),side=1,line=1.5,cex=pt9)
  if(i==1){mtext(expression(p[m]),side=2,line=1.5,cex=pt9);
      mtext("Seed addition",cex=pt10,adj=.5,outer=FALSE,side=2,line=2.5)}

  box()

  mtext(side=3, line =0,paste0(panelLabs[i]), cex = pt8 , adj = 0)

}
dev.off()



# - Legend ----

tiff(filename=paste0(primaryDirectory,"products/figures/identifiability-legend.tif"),
     height=3.5,width=6.5,units="in",res=300,compression="lzw",pointsize=12)

par(mfrow=c(1,1))
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend(x=0,y=0.5,c("Seed bag burial experiment", "Seed addition experiment"),
       pch=c(21,21),
       pt.bg=c("black","white"),
       box.lty=0, horiz=TRUE,
       text.width = c(1,.4), cex = pt10)
dev.off()
