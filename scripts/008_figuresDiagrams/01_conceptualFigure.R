####
####
# Script to create diagram of seed bag burial and addition experiments
####
####


# set font sizes
pt12 = 1
pt10 = 10/12
pt9 = 9/12
pt8 = 8/12
pt7 = 7/12
pt6 = 6/12

# - Seed bag burial experiment ----

tiff(filename="~/Dropbox/chapter-4/analysis/products/figures/seed-bag-trials.tif",
     height=1800,width=1600,units="px",res=800,compression="lzw",pointsize=12)

par(mfrow=c(1,1),mar=c(0,0,0,0),mgp=c(0,0.25,0),
    oma=c(2,1,0,0))
t.sample = t = c(0,12,24,36)/60

# - +Plot dimensions and axes ----
plot(NA,
     xlim=c(-1,40),ylim=c(0,30),
     type='n',frame=FALSE, cex = pt9,
     xlab = "Time (months)",
     ylab = "",
     axes=FALSE)

# - +Time axis ----
axis(1, c(0,10,20,30,40), col.ticks = 1,cex.axis = pt8,line = 0)
mtext("Time (months)",side=1,line=1,cex=pt9)

# - +Segments for trial length ----

arrows(0, 24, 12, 24, length=0.05, angle=90, code=1,lwd=1.25)
arrows(0, 21, 12, 21, length=0.05, angle=90, code=1,lwd=1.25)
arrows(0, 18, 12, 18, length=0.05, angle=90, code=1,lwd=1.25)

arrows(0, 15, 24, 15, length=0.05, angle=90, code=1,lwd=1.25)
arrows(0, 12, 24, 12, length=0.05, angle=90, code=1,lwd=1.25)
arrows(0, 9, 24, 9, length=0.05, angle=90, code=1,lwd=1.25)

arrows(0, 6, 36, 6, length=0.05, angle=90, code=1,lwd=1.25)
arrows(0, 3, 36, 3, length=0.05, angle=90, code=1,lwd=1.25)
arrows(0, 0, 36, 0, length=0.05, angle=90, code=1,lwd=1.25)

# - +Sampling points ----
points(x=rep(rev(c(12,24,36)),times=c(3,3,3)),
       y=seq(0,24,by=3),
       pch=c(21),
       bg=c("black"),
       col=c("black"))

# - +Time index ----
points(x=c(12,24,36),
       y=rep(26.5,3),
       pch=c(22),cex=2,
       bg=c("gray90"),
       col=c("gray90"))

text(x=c(12,24,36),y=rep(26.5,3),c(1:3),col='black',cex=pt7)
mtext("Time (index)",side=3,line=-1,cex=pt10)

# - +Bag index ----
rect(xleft=rep(-5,20),xright=rep(-.5,20),
     ybottom=seq(0,24,by=3)-1,ytop=seq(0,24,by=3)+1,
     border='gray90',col='gray90')

text(x=rep(-1.5,20),y=seq(0,24,by=3),
     rev(c(12,21,30,33,45,46,61,69,79)),
     col='black',cex=pt6)
mtext("Bag (index)",side=2,line=0,cex=pt9,outer=TRUE,adj = .4)

# - +'Legend' ----
segments(24,15,27,18.5,lty='dotted')
rect(xleft=16,xright=36,
     ybottom=18.5,ytop=23,lty='dotted')
text(x=26,y=21,
     "Intact seed and\n germinant counts",
     cex=pt7)

dev.off()


# - Seed addition experiment ----

tiff(filename="~/Dropbox/chapter-4/analysis/products/figures/seed-addition-trials.tif",
     height=1800,width=1600,units="px",res=800,compression="lzw",pointsize=12)

par(mfrow=c(1,1),mar=c(0,0,0,0),mgp=c(0,0.25,0),
    oma=c(2,1,0,0))
t.sample = t = c(0,12,24,36)/60

# - +Plot dimensions and axes ----
plot(NA,
     xlim=c(-1,40),ylim=c(0,7),
     type='n',frame=FALSE,
     xlab = "Time (months)",
     ylab = "",
     axes=FALSE)

# - +Time axis ----
axis(1, c(0,10,20,30,40), col.ticks = 1,cex.axis = pt8, line = 0)
mtext("Time (months)",side=1,line=1,cex=pt9)

# - +Segments for trial length----
arrows(0, 5, 36,  5, length=0.05, angle=90, code=1,lwd=1.25)
arrows(0, 3, 36,  3, length=0.05, angle=90, code=1,lwd=1.25)
arrows(0, 1, 36, 1, length=0.05, angle=90, code=1,lwd=1.25)

# - +Sampling points----
points(x=c(rep(c(12,24,36),3)),
       y=c(rep(c(1,3,5),each=3)),
       pch=c(21),
       bg=c("black"),
       col=c("black"))

# - +Time index----
points(x=c(12,24,36),
       y=rep(6,3),
       pch=c(22),cex=2,
       bg=c("gray90"),
       col=c("gray90"))

text(x=c(12,24,36),y=rep(6,3),c(1:3),col='black', cex = pt7)
mtext("Time (index)",side=3,line=-1,cex=pt10)

# - +Plot index----
rect(xleft=rep(-5,20),xright=rep(-.5,20),
     ybottom=seq(1,5,by=2)-.25,ytop=seq(1,5,by=2)+.25,
     border='gray90',col='gray90')

text(x=rep(-1.5,4),y=seq(1,5,by=2),
     c(7,19,21),
     col='black',cex=pt6)
mtext("Plot (index)",side=2,line=0,cex=pt9,adj=.4)


# - +'Legend'----
segments(36,5,27,4.5,lty='dotted')
rect(xleft=13,xright=32,
     ybottom=3.5,ytop=4.5,lty='dotted')
text(x=22,y=4,
     "Seedling counts",
     cex=pt7)


dev.off()
