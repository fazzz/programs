source("~/Rspa/ini.R")
source("~/Rspa/set_enecon.sh")
source("~/Rspa/set_res.R")
source("~/Rspa/plmGeneral_wsrange.R")

leg.pos="topright"
is.leg <- 0

xrange <- c(0,7,7)
yrange <- c(-7,0,7)

name.protein <- c("UBIQ","BPTI","SPE10")
num.protein <- length(name.protein)


dirbase <- "/home/yamamori/calspa/TAMD/"

for ( k in 1:4  ) {
  sen[k]<-c(3)
  n<-2+k
  iro[k]<-c(n)
  id.ys[k]<-c(2)
  tenshu[k]<-c(20)
}

par(mfrow=c(1,3))
par(mar=c(0.0,0.0,0.0,0.0))
par(oma=c(7.5,4.0,2.0,2.0))

for (i in 1:num.protein) {
  file.names[1]=paste(dirbase.tamd,name.protein[i],sep='')
  file.names[2]=paste(dirbase.tamd,name.protein[i],sep='')
  file.names[3]=paste(dirbase.tamd,name.protein[i],sep='')
  file.names[4]=paste(dirbase.tamd,name.protein[i],sep='')

  plmGeneralwsrange(data.names=file.names,
                    sd.names=file.names,
                    id.ys=id.ys,
                    ids.ys=ids.ys,
                    label.size=0.5,axis.size=2.0,
                    iro=iro,axis.ft="F",is.header="T",
                    sdiro=iro,
                    xrange=xrange,yrange=yrange,
                    sdyrange=yrange,
                    is.sen=sen,width=10.0,
                    ,warrow="T")
  text(3,-1,name.protein[i])
  box(lwd=2.0)
  if (i==1) {
    legend(3,-5.0,legend=c("vV-H'"),bty="n",
           lty=c(1),cex=1.0,y.intersp=2,
           col=iro[3],ncol=leg.col)
    legend(3,-4.0,legend=c("vV-H"),bty="n",
           lty=c(1),cex=1.0,y.intersp=2,
           col=iro[4],ncol=leg.col)
    legend(3,-3.0,legend=c("pc6-H"),bty="n",
           lty=c(1),cex=1.0,y.intersp=2,
           col=iro[1],ncol=leg.col)
    legend(3,-2.0,legend=c("pc6-H'"),bty="n",
           lty=c(1),cex=1.0,y.intersp=2,
           col=iro[2],ncol=leg.col)
  }

  if (i==1 ) {
    axis(2,yaxp=yrange,lwd=2.0)
  }

  axis(1,xaxp=xrange,lwd=2.0)
  mtext(outer=T,label.x,side=1,line=2.5,cex=0.8)
  if (i==1) {
    mtext(outer=T,label.y,side=2,line=2.5,cex=0.8)
  }
}
name.out=paste("~/papers/TAMD/fig3.econ_proteins",sep='')
OutTiff(name.out,w_val=880,h_val=330)
