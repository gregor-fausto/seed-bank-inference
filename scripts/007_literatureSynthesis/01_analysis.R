####
####
# Analyze literature synthesis for supplement
####
####

# - Libraries ----
library(tidyverse)

# - Read data ----
# table of literature survey
df <- read.csv("../data/01_literatureSynthesis/table-literature-survey.csv",header=TRUE)

# - Exploratory plots ----
# - + Seed bag data ----
df.sb  <-  df %>%
  dplyr::filter(experiment=="seed bag burial") %>%
  dplyr::mutate(years=totalMonths/12)

ggplot(df.sb) +
  geom_jitter(aes(x=seedNumber,y=repNumber),alpha=0.25) +
  theme_bw()

ggplot(df.sb) +
  geom_jitter(aes(x=repNumber,y=observationNumber),alpha=0.5) +
  theme_bw() +
  xlab("Number of replicate bags per observation") +
  ylab("Number of observations")

ggplot(df.sb) +
  geom_jitter(aes(x=repNumber,y=years),alpha=0.5) +
  theme_bw() +
  xlab("Number of replicate bags per observation") +
  ylab("Length of seed bag burial experiment (years)")

ggplot(df.sb) +
  geom_jitter(aes(x=years,y=observationNumber),alpha=0.5,height=0.1,width=0.1) +
  theme_bw()

hist(df.sb$observationNumber,breaks=10)

ggplot(df.sb) +
  geom_jitter(aes(x=years,y=observationNumber,size=seedNumber*repNumber),alpha=0.5) +
  theme_bw()

# - + Seed bag data and seed addition data ----
df2 <- df %>% dplyr::filter(experiment %in% c("seed addition","seed bag burial"))

ggplot(df2) +
  geom_jitter(aes(x=totalMonths,y=observationNumber),alpha=0.5) +
  facet_wrap(~experiment) +
  theme_bw()

ggplot(df2) +
  geom_jitter(aes(x=totalMonths,y=observationNumber,size=seedNumber*repNumber,color=experiment),alpha=0.5) +
  theme_bw()

# - Plots for supplement ----

df %>% 
  dplyr::filter(experiment %in% c("seed addition","seed bag burial")) %>%
  dplyr::group_by(experiment) %>%
  dplyr::summarise(n())

df <- df %>% 
  dplyr::filter(experiment %in% c("seed addition","seed bag burial")) 

dev.off()

# - +Summarize data ----

years <- c(1990,1995,2000,2005,2010,2015,2020)
summaryMat <- matrix(NA,nrow=6,ncol=2)
for(i in 1:(length(years)-1)){
  dfSub <- df[df$pubDate %in% (years[i]:years[i+1]),]
  dfVec <- dfSub %>% dplyr::group_by(experiment) %>%
    dplyr::summarise(n=n())
  summaryMat[i,] <- dfVec$n
}

dfSummary <- df %>% 
  dplyr::group_by(experiment,pubDate) %>%
  dplyr::summarise(n=n())

dfBag <- dfSummary[dfSummary$experiment=="seed bag burial",] %>% ungroup
dfAddition <- dfSummary[dfSummary$experiment=="seed addition",] %>% ungroup

# - ++ Output parameters ----
pt12 = 1
pt10 = 10/12
pt9 = 9/12
pt8 = 8/12
pt7 = 7/12

par(oma = c(2.5,4,0,0) + 0.1,
    mar = c(1,0,1,1) + 0.1)

# - +Panel A ----
# plots describing publication data

pdf("../products/figures/figure-meta-1.pdf",width=6,height=6)

par(mfrow=c(1,1))
plot(NA,xlim=c(1990,2021),
     ylim=c(0,6),
     xlab="Publication year",main='',ylab="Number of papers",cex=pt10,cex.axis=pt12,cex.lab=pt12)
rect(xleft=dfBag$pubDate-.45,xright=dfBag$pubDate,
     ybottom=0,ytop=dfBag$n,col='black',border='white')
rect(xleft=dfAddition$pubDate,xright=dfAddition$pubDate+.45,
     ybottom=0,ytop=dfAddition$n,col='orange',border='white')
legend("topleft",bty='n',cex=pt10,
       c("Demographic model parameterized with seed bag burial experiments",
       "Demographic model parameterized with seed addition experiments"),
       pch=15,col=c('black','orange'))

dev.off()

# - +Summarize # observations ----
dfSummary=df %>% 
  dplyr::group_by(experiment,observationNumber) %>%
  dplyr::summarise(n=n())

dfBag = dfSummary[dfSummary$experiment=="seed bag burial",] %>% ungroup
dfAddition = dfSummary[dfSummary$experiment=="seed addition",] %>% ungroup

pdf("../products/figures/figure-meta-2.pdf",width=6,height=6)

# - +Panel B ----
#  plots describing publication data
par(mfrow=c(1,1))
plot(NA,
     ylim=c(0,14),xlim=c(0,20),
     xlab="Number of experimental observations",main='',ylab="Number of papers",cex=pt10,cex.axis=pt12,cex.lab=pt12)
rect(xleft=dfBag$observationNumber-.45,xright=dfBag$observationNumber,
     ybottom=0,ytop=dfBag$n,col='black',border='white')
rect(xleft=dfAddition$observationNumber,xright=dfAddition$observationNumber+.45,
     ybottom=0,ytop=dfAddition$n,col='orange',border='white')

dev.off()

# - +Split by experiment ----
dfBag = df[df$experiment=="seed bag burial",]
dfAddition = df[df$experiment=="seed addition",]

# - +Panel C ----
#  summarize experiments by length vs. # observations

pdf("../products/figures/figure-meta-3.pdf",width=6,height=6)

par(mfcol=c(1,1))
plot(dfBag$totalMonths,dfBag$observationNumber,
     pch=16,xlim=c(0,140),ylim=c(0,20),type='n',
     xlab="Length of experiment (months)",
     ylab="Number of experimental observations",cex=1,cex.axis=pt12,cex.lab=pt12)
abline(a=0,b=1/12,lty='dashed')
points(dfAddition$totalMonths,dfAddition$observationNumber,pch=21,bg='white',col='orange',cex=1.5)
points(dfBag$totalMonths,dfBag$observationNumber,pch=21,col='black',cex=1)

dev.off()

# - +Panel D ----
#  summarize by sample size vs # observations
pdf("../products/figures/figure-meta-4.pdf",width=6,height=6)

par(mfcol=c(1,1))
plot(dfBag$repNumber*dfBag$seedNumber,dfBag$observationNumber,
     pch=16,type='n',xlim=c(0,9000),ylim=c(0,20),
     xlab="Total number of seeds per observation (bags times seeds/bag)",
     ylab="Number of experimental observations",cex.axis=pt12,cex.lab=pt12)
points(dfAddition$repNumber*dfAddition$seedNumber,dfAddition$observationNumber,pch=21,bg='white',col='orange',cex=1.5)
points(dfBag$repNumber*dfBag$seedNumber,dfBag$observationNumber,pch=21,col='black',cex=1)

dev.off()