#!/bin/zsh

DataName=Data.txt
Xid=2;Yid=3
OutName=Tst
cat << eof > dihfel.R
Data <- read.table("$DataName",header=T)
label.x <- expression(paste(theta[1],'(deg.)',sep=''))
label.y <- expression(paste(theta[2],'(deg.)',sep=''))
NumFrm <- length(Data[,1])
Bin <- 1
NumBin <- 360 / Bin
Density <- array(rep(0,NumBin*NumBin),c(NumBin,NumBin))
Dih <- array(seq(1:2*NumBin),c(2,NumBin))
for (i in 1:NumBin) {
  Dih[1,i] <- -180 - Bin / 2 + Bin * i
  Dih[2,i] <- -180 - Bin / 2 + Bin * i
}
Xid <- $Xid
Yid <- $Yid
for (i in 1:NumFrm) {
  i1 <- (Data[i,Xid]+180) / Bin + 1  
  i2 <- (Data[i,Yid]+180) / Bin + 1  
  Density[i1,i2] <- Density[i1,i2] + 1
}
MaxVal <- 0
for (i in 1:NumBin) {
  for (j in 1:NumBin) {
    if (Density[i,j] > 0 ) Density[i,j] <- -log(Density[i,j]/NumFrm)
    if (Density[i,j] > MaxVal ) MaxVal <- Density[i,j]
  }
}
MinVal <- 10**5
for (i in 1:NumBin) {
  for (j in 1:NumBin) {
    if (Density[i,j] == 0 ) Density[i,j] <- MaxVal
    if (Density[i,j] < MinVal ) MinVal <- Density[i,j]
  }
}
cat(MaxVal, MinVal,'\n')
for (i in 1:NumBin) {
  for (j in 1:NumBin) {
    Density[i,j] <- Density[i,j] - MinVal
  }
}
MyColor <- function(n,alpha=1)
{
    if ((n <- as.integer(n[1L])) > 0) {
        j <- n%/%3
        k <- n%/%3
        i <- n - j - k
        c(if (i > 0) hsv(h = seq.int(from = 40/60, to = 25/60, 
            length.out = i), alpha = alpha), if (j > 0) hsv(h = seq.int(from = 23/60, 
            to = 11/60, length.out = j), alpha = alpha), if (k > 
            0) hsv(h = seq.int(from = 8/60, to = 0/60, length.out = k-1), 
            alpha = alpha, s = seq.int(from = 1, to = 0.9, length.out = k-1), 
            v = 1),hsv(0,0,1))
    }
    else character(0L)
}
file.name <- '${OutName}.tiff'
tiff(file.name,width=1000,height=1000)
#par(mfrow = c(1,4))
par(oma = c(0,0,0,0) )
source("FillConMap.R")
source("FillConBar.R")
xrange <- c(-180,180,6)
yrange <- c(-180,180,6)
nf <- layout(matrix(c(1,2,4,3,3,4),2,3,byrow=TRUE),c(4,4,1),c(4,5))
par(mar =c(5,5,2,2) )
FillConMap(Dih[1,],Dih[2,],xlab=label.x,ylab=label.y,Density,xrange[1:2],yrange[1:2],color=MyColor,nlevels=10)
FillConMap(Dih[1,],Dih[2,],xlab=label.x,ylab=label.y,Density,xrange[1:2],yrange[1:2],color=MyColor,nlevels=6)
FillConMap(Dih[1,],Dih[2,],xlab=label.x,ylab=label.y,Density,xrange[1:2],yrange[1:2],color=MyColor,nlevels=4)
FillConBar(Dih[1,],Dih[2,],xlab=label.x,ylab=label.y,Density,xrange[1:2],yrange[1:2],color=MyColor,nlevels=10)
eof
Rscript dihfel.R
display $OutName.tiff
