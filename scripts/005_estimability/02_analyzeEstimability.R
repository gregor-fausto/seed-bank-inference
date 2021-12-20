####
####
# Study parameter estimability in non-parametric models
# this is for a constant mortality constant germination model
# with different number of time points
####
####

# - Parameters ----
allParams = c("p.m","p.g")

for(k in 1:2){
  
  # remove all but the index and the vector of parameters
  rm(list=setdiff(ls(all=TRUE), c("k","allParams"))) 
  
  focalParam = allParams[k]
  
  
  
  # - Libraries ----
  library(MCMCvis)
  library(tidyverse)
  source("scripts/005_estimability/01_functionsEstimability.R")
  
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
  
  # get names of directories
  tempDirs=list.dirs("outputs/002_statisticalModelFitting/02_estimability",recursive=FALSE)
  
  # construct parameter table from file names
  parm.table=data.frame(do.call(rbind,lapply(tempDirs,f)))
  
  # - +Seed bag burial experiment ----
  
  # list to hold evaluation output
  eval.list = list()
  
  # - For loop to evaluate models ----
  
  for(i in 1:dim(parm.table)[1]){
    
    # - +Extract posteriors ----
    
    # get current directory
    tmp = tempDirs[i]
    parm.table.tmp = parm.table[i,]
    
    # get filenames in directories
    fileNamesBags=list.files(paste0(tmp,"/posteriors-NpCmCg/seedBagBurial/"))
    fileNamesAddition=list.files(paste0(tmp,"/posteriors-NpCmCg/seedAddition/"))
    
    # list of all file paths for all posteriors
    fitsBags = paste0(paste0(tmp,"/posteriors-NpCmCg/seedBagBurial/",fileNamesBags))
    fitsAddition = paste0(paste0(tmp,"/posteriors-NpCmCg/seedAddition/",fileNamesAddition))
    
    # if directory is empty (no files), go to next one
    if(length(fitsBags)<1) next
    
    # extract posteriors for mortality and germination
    posterior.bags.all = lapply(fitsBags,f.extract,params=c("p.m","p.g"))
    posterior.addition.all = lapply(fitsAddition,f.extract,params=c("p.m","p.g"))
    
    # - +Model assessment ----
    
    eval.bags = wrapper.fun(posterior.list = posterior.bags.all, par = focalParam, parm.df = parm.table.tmp)
    eval.addition = wrapper.fun(posterior.list = posterior.addition.all, par = focalParam, parm.df = parm.table.tmp)
    
    # add experiment variable
    eval.out = rbind(data.frame(c(eval.bags,experiment="bags")),data.frame(c(eval.addition,experiment="addition")))
    
    # put in list
    eval.list[[i]] = eval.out
    
  }
  
  # - Bind data frame ----
  
  df = do.call("rbind", eval.list)
  
  
  # - Plotting ----
  
  panels <- unique(df[,c("p.m","p.g")])
  
  # set font sizes
  pt12 = 1
  pt10 = 10/12
  pt9 = 9/12
  pt8 = 8/12
  pt7 = 7/12
  
  # - +Panels for main text ----
  
  # - ++Bias ----
  
  tiff(filename=paste0("products/figures/summary-",focalParam,".tif"),
       height=3.5,width=6,units="in",res=300,compression="lzw",pointsize=12)
  
  if(focalParam=="p.m"){panelLabs =  c("A.","B.","C.","D.")} else {panelLabs = c("I.", "J.", "K.", "L.")}
  
  par(mfrow = c(2,4),
      oma = c(2.5,4,1.1,0) + 0.1,
      mar = c(0.4,0,1,.4) + 0.1)
  
  for(i in 1:dim(panels)[1]){
    
    # get data for panel
    df.tmp <- df %>% dplyr::filter(p.m==panels$p.m[i] & p.g==panels$p.g[i])
    
    if(focalParam=="p.m"){yax = .02}else{yax = 0}
    plot(NA,
         xlim=c(0,32), ylim = c(-.02, .02+yax),
         ylab='',
         xlab="",
         yaxt='n',
         xaxt='n')
    axis(1, at = seq(0 , 30, by = 5), labels=F)
    
    if(i==1){ axis(2, at = seq(-.02, .02+yax, by = .01), labels = F);
      axis(2, at = seq(-.02, .02+yax, by = .02), las=1, cex = pt7) }
    
    abline(h = 0, col = "gray50")
    box()
    
    segments(x0=df.tmp$n.bags + c(-0.75,0.75),
             y0=df.tmp$bias.mean-df.tmp$bias.ci95,
             y1=df.tmp$bias.mean+df.tmp$bias.ci95)
    points(x= df.tmp$n.bags+ c(-0.75,0.75), y = df.tmp$bias.mean ,
           pch=ifelse(df.tmp$experiment=="bags",21,21),
           bg=ifelse(df.tmp$experiment=="bags","black","white"))
    
    mtext(side=3, line =0,paste0(panelLabs[i]), cex = pt8 , adj = 0)
    if(i==1){ mtext(side=3, line =0,expression(paste(p[m]==0.1, ",  ",  p[g]==0.1)), cex = pt8 , adj = 0.5) }
    if(i==2){ mtext(side=3, line =0,expression(paste(p[m]==0.1, ",  ",  p[g]==0.5)), cex = pt8 , adj = 0.5) }
    if(i==3){ mtext(side=3, line =0,expression(paste(p[m]==0.5, ",  ",  p[g]==0.1)), cex = pt8 , adj = 0.5) }
    if(i==4){ mtext(side=3, line =0,expression(paste(p[m]==0.5, ",  ",  p[g]==0.5)), cex = pt8 , adj = 0.5) }
    
    
    if(i==1){
      mtext(side=2, line = 3, "Bias", cex = pt9)
    }
    
  }
  
  # - ++Width ----
  
  if(focalParam=="p.m"){panelLabs =  c("E.","F.","G.","H.")} else {panelLabs = c("M.", "N.", "O.", "P.")}
  
  for(i in 1:dim(panels)[1]){
    
    # get data for panel
    df.tmp <- df %>% dplyr::filter(p.m==panels$p.m[i] & p.g==panels$p.g[i])
    
    if(focalParam=="p.m"){yax = .2}else{yax = 0}
    plot(NA,
         xlim=c(0,32), ylim = c(0, .2+yax),
         ylab='',
         xlab="",
         yaxt='n',
         xaxt='n')
    axis(1, at = seq(0 , 30, by = 10), las=1, cex = pt7,padj=-.5)
    if(i==1){  axis(2, at = seq(0, .2+yax, by = .05), labels=F );
      axis(2, at = seq(0, .2+yax, by = .1), las=1, cex = pt7) }
    
    abline(h = 0, col = "gray50")
    box()
    
    segments(x0=df.tmp$n.bags + c(0,0),
             y0=df.tmp$width.mean-df.tmp$width.ci95,
             y1=df.tmp$width.mean+df.tmp$width.ci95)
    points(x= df.tmp$n.bags+ c(0,0), y = df.tmp$width.mean ,
           pch=ifelse(df.tmp$experiment=="bags",21,21),
           bg=ifelse(df.tmp$experiment=="bags","black","white"))
    
    mtext(side=3, line =0,paste0(panelLabs[i]), cex = pt8 , adj = 0)
    
    if(i==1){
      mtext(side=2, line = 3, "95% credible interval width", cex = pt9)
    }
    
  }
  mtext(side=1, line = 1.25, "Number of samples (bags or plots)", cex = pt9,  outer = TRUE)
  if(focalParam=="p.m"){
    mtext(side=3, line =0, expression(Probability~of~mortality~p[m]), adj=-.125, cex = pt10,  outer = TRUE)
  } else {
    mtext(side=3, line =-.1, expression(Probability~of~germination~p[g]), adj=-.125, cex = pt10,  outer = TRUE)
  }
  
  
  dev.off()
  
  
  # - +Panels for supplement ----
  # - ++Coverage ----
  
  tiff(filename=paste0("products/figures/summary-",focalParam,"-supplement.tif"),
       height=3.5,width=6,units="in",res=300,compression="lzw",pointsize=12)
  
  panelLabs =  c("A.","B.","C.","D.","I.", "J.", "K.", "L.")
  
  par(mfrow = c(2,4),
      oma = c(2.5,4,1.1,0) + 0.1,
      mar = c(0.4,0,1,.4) + 0.1)
  
  for(i in 1:dim(panels)[1]){
    
    # get data for panel
    df.tmp <- df %>% dplyr::filter(p.m==panels$p.m[i] & p.g==panels$p.g[i])
    
    plot(NA,
         xlim=c(0,32), ylim = c(.75, 1),
         ylab='',
         xlab="",
         yaxt='n',
         xaxt='n')
    axis(1, at = seq(0 , 30, by = 10), labels=F)
    if(i==1){ axis(2, at = seq(.75, 1, by = .025), labels = F);
      axis(2, at = seq(.75, 1, by = .05), las=1, cex = pt7) }
    
    abline(h = .95, col='gray50')
    box()
    
    segments(x0=df.tmp$n.bags+ c(-0.75,0.75),
             y0=df.tmp$coverage.lo,
             y1=df.tmp$coverage.hi)
    points(x= df.tmp$n.bags+ c(-0.75,0.75), y = df.tmp$coverage.mean ,
           pch=ifelse(df.tmp$experiment=="bags",21,21),
           bg=ifelse(df.tmp$experiment=="bags","black","white"))
    
    mtext(side=3, line =0,paste0(panelLabs[(1:4)[i]]), cex = pt8 , adj = 0)
    if(i==1){ mtext(side=3, line =0,expression(paste(p[m]==0.1, ",  ",  p[g]==0.1)), cex = pt8 , adj = 0.5) }
    if(i==2){ mtext(side=3, line =0,expression(paste(p[m]==0.1, ",  ",  p[g]==0.5)), cex = pt8 , adj = 0.5) }
    if(i==3){ mtext(side=3, line =0,expression(paste(p[m]==0.5, ",  ",  p[g]==0.1)), cex = pt8 , adj = 0.5) }
    if(i==4){ mtext(side=3, line =0,expression(paste(p[m]==0.5, ",  ",  p[g]==0.5)), cex = pt8 , adj = 0.5) }
    
    if(i==1){
      mtext(side=2, line = 3.1, "Coverage", cex = pt9)
    }
    
    mtext(side=1, line = 2.5, "Number of samples (bags or plots)", cex = pt9,  outer = TRUE)
    
  }
  
  # - ++RMSE ----
  
  for(i in 1:dim(panels)[1]){
    
    # get data for panel
    df.tmp <- df %>% dplyr::filter(p.m==panels$p.m[i] & p.g==panels$p.g[i])
    
    if(focalParam=="p.m"){yax = .01}else{yax = 0}
    plot(NA,
         xlim=c(0,32), ylim = c(0, .01+yax),
         ylab='',
         xlab="",
         yaxt='n',
         xaxt='n')
    axis(1, at = seq(0 , 30, by = 10), las=1, cex = pt7,padj=-.5)
    if(i==1){  axis(2, at = seq(0, .01+yax, by = .0025), labels=F );
      axis(2, at = seq(0, .01+yax, by = .005), las=1, cex = pt7) }
    
    abline(h = 0, col = "gray50")
    box()
    
    points(x= df.tmp$n.bags+ c(-0.5,0.5), y = df.tmp$rmse ,
           pch=ifelse(df.tmp$experiment=="bags",21,21),
           bg=ifelse(df.tmp$experiment=="bags","black","white"))
    
    if(i==1){
      mtext(side=2, line = 3.1, "Root mean squared error", cex = pt9)
    }
    
  }
  
  mtext(side=3, line =0,paste0(panelLabs[(5:8)[i]]), cex = pt8 , adj = 0)
  mtext(side=1, line = 1.25, "Number of samples (bags or plots)", cex = pt9,  outer = TRUE)
  if(focalParam=="p.m"){
    mtext(side=3, line =0, expression(Probability~of~mortality~p[m]), adj=-.125, cex = pt10,  outer = TRUE)
  } else {
    mtext(side=3, line =-.1, expression(Probability~of~germination~p[g]), adj=-.125, cex = pt10,  outer = TRUE)
  }
  
  dev.off()
  
  # - Legend ----
  
  tiff(filename="products/figures/summary-legend.tif",
       height=3.5,width=6.5,units="in",res=300,compression="lzw",pointsize=12)
  
  par(mfrow=c(1,1))
  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  legend(x=0,y=0.5,c("Seed bag burial experiment", "Seed addition experiment"),
         pch=c(21,21),
         pt.bg=c("black","white"),
         lty=c(1,1),box.lty=0, horiz=TRUE, 
         text.width = c(1,.4), cex = pt10)
  dev.off()
}