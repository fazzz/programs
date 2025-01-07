source("~/Rspa/ini.R")
source("~/Rspa/set_spe.R")
source("~/Rspa/set_res.R")
source("~/Rspa/plmGeneral_wsrange.R")

leg.pos="topright"
is.leg <- 0

label.x <- expression(paste("frequency / cm"^-1))

xrange <- c(0,4000,4)
#yrange <- c(0,0.2,2)

#name.protein <- c("UBIQ","BPTI","HP35")
name.protein <- c("UBIQ")
num.protein <- length(name.protein)

tp <- c("300")
num.tp <- length(tp)

dirbase <- "/home/yamamori/calspa/refcalc/"

sen=c(1,2)
iro<-c(3,2)
id.ys<-c(2,2)
tenshu<-c(1,1)

par(mfrow=c(1,1))
par(mar=c(0.0,0.0,0.0,0.0))
par(oma=c(7.5,4.0,2.0,2.0))

file.names<-NULL
for (i in 1:num.protein) {
  file.names[1]=paste(dirbase,name.protein[i],"/spevac/spectrum/spe_",name.protein[i],"v_300_1_c_a.txt",sep='')
  file.names[2]="/home/yamamori/calspa/refcalc/UBIQ/nmadas/UBIQ.freq_hist_40_4"

  if (i==1 ) {
    yrange <- c(0.4,30,1)
  }
  if (i==2 ) {
    yrange <- c(0.65,20,1)
  }
  if (i==3 ) {
    yrange <- c(0.6,18,1)
  }
  plmGeneralwsrange(data.names=file.names,
                    id.ys=id.ys,
                    label.size=0.5,axis.size=2.0,
                    iro=iro,axis.ft="F",
                    xrange=xrange,yrange=yrange,
                    is.sen=sen,width=40.0,is.header="F")
#  text(800,25,name.protein[i])
  box(lwd=2.0)

  if (i==1) {
    legend(700,0.1,legend=c("pc6"),bty="n",
           lty=c(1),cex=1.0,y.intersp=2,
           col=iro[1],ncol=leg.col)
    legend(700,0.13,legend=c("vV"),bty="n",
           lty=c(1),cex=1.0,y.intersp=2,
           col=iro[2],ncol=leg.col)
  }

  axis(1,xaxp=xrange,lwd=2.0,cex.axis=1.5)
  mtext(outer=T,label.x,side=1,line=5.0,cex=2.0)
#  if (i==1 ) {
#    axis(2,yaxp=yrange,lwd=2.0)
#  }
  if (i==1) {
    mtext(outer=T,label.y,side=2,line=2.5,cex=2.0)
  }
}

name.out=paste("~/papers/TAMD/fig1.biophys_spec_proteins_ref",sep='')
OutTiff(name.out,w_val=800,h_val=600)
