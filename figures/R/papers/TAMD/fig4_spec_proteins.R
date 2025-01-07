source("~/Rspa/ini.R")
source("~/Rspa/set_spe.R")
source("~/Rspa/set_res.R")
source("~/Rspa/plmGeneral_wsrange.R")

leg.pos="topright"
is.leg <- 0

label.x <- expression(paste("frequency / cm"^-1))
label.y <- paste("I(",expression(omega),")/N /cm")

xrange <- c(0,1000,10)
yrange <- c(0,0.2,2)

name.protein <- c("UBIQ","BPTI","SPE10")
num.protein <- length(name.protein)

tp <- c("300")
num.tp <- length(tp)

dirbase <- "/home/yamamori/calspa/TAMD/"

for ( k in 1:4  ) {
  sen[k]<-c(3)
  n<-2+k
  iro[k]<-c(n)
  id.ys[k]<-c(2)
  tenshu[k]<-c(20)
}

par(mfrow=c(3,3))
par(mar=c(0.0,0.0,0.0,0.0))
par(oma=c(7.5,4.0,2.0,2.0))

for (i in 1:num.protein) {
  file.names[1]=paste(dirbase,name.protein[i],"dummy",sep='')
  file.names[2]=paste(dirbase,name.protein[i],"dummy",sep='')

  plmGeneralwsrange(data.names=file.names,
                    id.ys=id.ys,
                    label.size=0.5,axis.size=2.0,
                    iro=iro,axis.ft="F",
                    xrange=xrange,yrange=yrange,
                    is.sen=sen,width=10.0,is.header="F")
  text(800,0.03,name.protein[i])
  box(lwd=2.0)

  if (i==1) {
    legend(700,0.1,legend=c("pc6"),bty="n",
           lty=c(1),cex=1.0,y.intersp=2,
           col=iro[1],ncol=leg.col)
    legend(700,0.13,legend=c("vV"),bty="n",
           lty=c(1),cex=1.0,y.intersp=2,
           col=iro[2],ncol=leg.col)
  }

  if (i==1 ) {
    axis(2,yaxp=yrange,lwd=2.0)
  }
  if (i==1) {
    mtext(outer=T,label.y,side=2,line=2.5,cex=0.8)
  }
}

for (i in 1:num.protein) {
  file.names[1]=paste(dirbase,name.protein[i],"dummy",sep='')
  file.names[2]=paste(dirbase,name.protein[i],"dummy",sep='')
  file.names[3]=paste(dirbase,name.protein[i],"dummy",sep='')

  plmGeneralwsrange(data.names=file.names,
                    label.names=label.names,
                    id.ys=id.ys,
                    label.size=0.5,axis.size=2.0,
                    iro=iro,axis.ft="F",
                    xrange=xrange,yrange=yrange,
                    is.sen=sen,width=10.0,is.header="F")
  text(800,0.06,name.protein[i])
  text(800,0.03,"pc6")
  box(lwd=2.0)

  if (i==1) {
    legend(650,0.1,legend=c("dt=1fs"),bty="n",
           lty=c(1),cex=1.0,y.intersp=2,
           col=iro[1],ncol=leg.col)
    legend(650,0.13,legend=c("dt=4fs"),bty="n",
           lty=c(1),cex=1.0,y.intersp=2,
           col=iro[2],ncol=leg.col)
    legend(650,0.16,legend=c("dt=5fs"),bty="n",
           lty=c(1),cex=1.0,y.intersp=2,
           col=iro[3],ncol=leg.col)
  }

  
  if (i==1 ) {
    axis(2,yaxp=yrange,lwd=2.0)
  }
}

for (i in 1:num.protein) {
  file.names[1]=paste(dirbase,name.protein[i],"dummy",sep='')
  file.names[2]=paste(dirbase,name.protein[i],"dummy",sep='')
  file.names[3]=paste(dirbase,name.protein[i],"dummy",sep='')

  if (i==1)
    label.names<-c('dt=1fs','dt=4fs','dt=5fs')
  else
    is.leg <- 0
  plmGeneralwsrange(data.names=file.names,
                    label.names=label.names,
                    id.ys=id.ys,
                    label.size=0.5,axis.size=2.0,
                    iro=iro,axis.ft="F",
                    xrange=xrange,yrange=yrange,
                    is.sen=sen,width=10.0,is.header="F")
  text(800,0.06,name.protein[i])
  text(800,0.03,"vV")
  box(lwd=2.0)

  if (i==1) {
    legend(650,0.1,legend=c("dt=1fs"),bty="n",
           lty=c(1),cex=1.0,y.intersp=2,
           col=iro[1],ncol=leg.col)
    legend(650,0.13,legend=c("dt=4fs"),bty="n",
           lty=c(1),cex=1.0,y.intersp=2,
           col=iro[2],ncol=leg.col)
    legend(650,0.16,legend=c("dt=5fs"),bty="n",
           lty=c(1),cex=1.0,y.intersp=2,
           col=iro[3],ncol=leg.col)
    legend(650,0.19,legend=c("dt=6fs"),bty="n",
           lty=c(1),cex=1.0,y.intersp=2,
           col=iro[4],ncol=leg.col)
  }

  if (i==1 ) {
    axis(2,yaxp=yrange,lwd=2.0)
  }
  axis(1,xaxp=xrange,lwd=2.0)
  mtext(outer=T,label.x,side=1,line=2.5,cex=0.8)
}

name.out=paste("~/papers/TAMD/fig4.spec_proteins",sep='')
OutTiff(name.out,w_val=880,h_val=660)
