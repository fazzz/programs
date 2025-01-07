source("~/Rspa/ini.R")
source("~/Rspa/set_spe.R")
source("~/Rspa/set_res.R")
source("~/Rspa/plmGeneral_wsrange.R")

leg.pos="topright"
is.leg <- 0

xrange <- c(0,500,5)
yrange <- c(0,0.9,1)

name.res <- c("GLY","GLN","GLU","HIE","ILE","PRO")

num.res <- length(name.res)

dirbase.name <- "/home/yamamori/calspa/TAMD/dipeptides/"

for ( k in 1:2  ) {
  n<-1+k
  id.ys[k]<-c(2)
  tenshu[k]<-c(20)
}

iro<-c(2,4)
sen<-c(1,1)
senshu<-c(1,2)
hutosa<-c(2.0,2.0)

par(mfrow=c(2,3))
par(mar=c(0.0,0.0,0.0,0.0))
par(oma=c(7.5,4.0,2.0,2.0))

file.names<-NULL
for (i in 1:num.res) {
  file.names[1]=paste(dirbase.name,name.res[i],"/spe_wdi_NVE_amber/spectrum/spe_",name.res[i],"Dv_",temp[1],"_1_ave_c_a.txt",sep='')
  file.names[2]=paste(dirbase.name,name.res[i],"/spe_wdi_NVE_amber/spectrum/spe_",name.res[i],"Dv_",temp[1],"_5_ave_c_a.txt",sep='')

  plmGeneralwsrange(data.names=file.names,
                      sd.names=file.names,
#                      label.names=label.names,
                      id.ys=id.ys,
                      ids.ys=ids.ys,
                      label.size=0.5,axis.size=2.0,
                      iro=iro,axis.ft="F",is.header="T",
                      sdiro=iro,
                      xrange=xrange,yrange=yrange,
                      sdyrange=yrange,
                      is.sen=sen,width=10.0,
                      ,warrow="F")
#  text(360,0.6,name.res[i],cex=2.0)
  box(lwd=2.0)

  if (i==1 ) {
    legend(400,1.0,legend=c("1fs"),bty="n",
           lty=c(1),cex=1.0,y.intersp=2,
           col=iro[1],ncol=leg.col)
    legend(400,1.0,legend=c("5fs"),bty="n",
           lty=c(1),cex=1.0,y.intersp=2,
           col=iro[2],ncol=leg.col)
  }
  
  if ( i==5 ) {
#    mtext(outer=T,label.x,side=1,line=5.0,cex=1.5)
  }
  if (i==4 || i==5 || i==6 ) {
    axis(1,xaxp=c(0,400,4),lwd=2.0,cex.axis=1.7)
  }
  
  if (i==1 || i==4  ) {
    axis(2,yaxp=c(0,0.5,1),lwd=2.0,cex.axis=1.7)
#    mtext(outer=T,label.y,side=2,line=2.5,cex=1.5)
  }

}
name.out=paste("~/papers/TAMD/fig.tyuukan_3_spe_con_6_aminoacids",sep='')
OutTiff(name.out,w_val=880,h_val=660)
