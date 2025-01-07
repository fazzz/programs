#!/bin/zsh
OutName=Fig2TimeRmsdComd
if [ -e rmsd.fl ];then rm rmsd.fl;fi
if [ -e comd.fl ];then rm comd.fl;fi
for i in Exp Atd Gld Mvd;do
  \ls /home/takemura/LyzTriNag/Analysis/Rmsd/TimeRmsd/Rmsd${i}CmpS1ns.txt >> rmsd.fl
  \ls /home/takemura/LyzTriNag/Analysis/ComDist/TimeComDist/ComDist${i}CmpS1ns.txt >> comd.fl
done
Rdir=/home/takemura/Rscript
RefData=/home/takemura/LyzTriNag/Analysis/ComDist/RefData/ComDistRef.txt
cat << eof > 1.R
file.name <- '${OutName}.eps'
postscript(file.name,width=3.5,height=5.2,horizontal=FALSE,onefile=FALSE,paper="special")
par(mfrow = c(2,1))
par(oma = c(3.5,4.5,0.5,0.5) )
par(mar =c(0,0,0,0) )
#par(mar =c(5,7,2,2) )
source('$Rdir/PlAxis.R')
source('$Rdir/SrchRange.R')
pttrn  <- c('brown',2:6,'gray50','orange','purple','indianred')
iro    <- rep(pttrn,100)
NumSys <- 4
SysNames <- c('Crystal','AutoDock','GOLD','MVD')
RmsdFiles <- readLines('rmsd.fl')
ComdFiles <- readLines('comd.fl')
# Data
RefData <- read.table('$RefData',header=T)
RefComDist <- RefData[,2]
a <- read.table(RmsdFiles[1],header=T)
Ns <- a[,1]
NumFrm <- length(Ns)
Rmsd <- array(seq(1:(NumFrm*NumSys)),c(NumSys,NumFrm))
Comd <- array(seq(1:(NumFrm*NumSys)),c(NumSys,NumFrm))
for (i in 1:NumSys) {
  a <- read.table(RmsdFiles[i],header=T)
  Rmsd[i,] <- a[,2]
  a <- read.table(ComdFiles[i],header=T)
  Comd[i,] <- a[,2]
}
label.x <- expression(paste('Time (ns)',sep=''))
xrange <- SrchRange(Ns,5)
Evy <- 2
# (1) RMSD
label.y <- expression(paste(RMSD[Complex],' (',ring(A),')',sep=''))
yrange <- SrchRange(c(min(Rmsd),max(Rmsd)),0.5)
PlYAxis(1.0)
for (i in 1:NumSys) {
  lines(Ns[seq(Evy,NumFrm,Evy)],Rmsd[i,][seq(Evy,NumFrm,Evy)],col=iro[i])
  x.pos <- xrange[1] + (xrange[2] - xrange[1]) * 0.2 * i
  y.pos <- yrange[1] + (yrange[2] - yrange[1]) * 0.1
  text(x.pos,y.pos,SysNames[i],col=iro[i],cex=1.0)
}
mtext(label.y,2,3,cex=1.2)
# (2) Comd
label.y <- expression(paste(Delta,R[COM],' (',ring(A),')',sep=''))
#yrange <- SrchRange(c(min(Comd),max(Comd)),2)
yrange <- c(9,17,4)
PlAxis(1.0)
for (i in 1:NumSys) {
  lines(Ns[seq(Evy,NumFrm,Evy)],Comd[i,][seq(Evy,NumFrm,Evy)],col=iro[i])
}
par(xpd=TRUE)
for (i in 1:NumSys) {
#  lines(xrange[1:2],c(RefComDist[i],RefComDist[i]),lty=5,lw=3,col=iro[i])
  lines(c(-1,21),c(RefComDist[i],RefComDist[i]),lty=5,lw=3,col=iro[i])
}
mtext(label.y,2,2.8,cex=1.2)
mtext(label.x,1,1.8,cex=1.2)
eof
Rscript 1.R
display ${OutName}.eps
