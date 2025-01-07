setwd("~/Rspa")
source("~/Rspa/ini.R")
source("~/Rspa/plmGeneral_wsrange_wdrange_wdiffx.R")

proname<-"AD"

TLbase<-"10"
numEX<-"250"
numREUS<-"16"
pname<-"UmbAt_Nbin_4x4_K_1"

name.title <- paste("/home/yamamori/gakkai/BioPhys_2013/eps/fig_REUS_AD_dihedtraj_2013-10-20",sep="")

leg.pos="topright"
is.leg <- 0

fact.y=1/pi

file.name <- paste(name.title,".eps",sep="")
postscript(file.name,width=1.719,height=2.13444,horizontal=FALSE,onefile=FALSE,paper="special")

par(mfrow=c(2,1))
par(mar=c(0.0,0.0,0.0,0.0))
par(oma=c(0.5,0.5,0.1,0.1))

par(cex=1.0)

filex.names<-NULL
filey.names<-NULL

xrange <- c(0,10000,5)
xrange.axis <- c(0,10000,5)

yrange <- c(-1.0,1.0,4)
yrange.axis <- c(-1.0,1.0,4)

label.x <- expression(paste("time (ps)"))
label.y <- expression(paste("Dihedral Angle"))

dir <- "/home/yamamori/calspa/refcalc/REUS/AD/s_REUSVAC_2013-08-12_ff99SB/"

hutosa <- NULL

for (i in 1:2) {
    sen<-c(0)
    senshu <- 1
    iro<-c(2)
    id.ys<-i+1

    dirbase.name <- paste(dir,"/",pname,"/nEX_",numEX,"/freq_",TLbase,"ps/",sep="")
        
    filex.names<-paste(dirbase.name,"anl/","RST_vs_t.1_1_2_n",sep='')
    filey.names<-paste(dirbase.name,"anl/","RST_vs_t.1_1_2_n",sep='')

    hutosa <- 1.0
    ave.x <- 1
          
    plmGeneralwsrangewdrangewdiffx(datax.names=filex.names,
                                   datay.names=filey.names,
                                   sd.names=file.names,
                                   id.ys=id.ys,
                                   ids.ys=ids.ys,
                                   label.size=0.5,axis.size=2.0,
                                   iro=iro,axis.ft="F",is.header="F",
                                   sdiro=iro,
                                   xrange=xrange,yrange=yrange,
                                   sdyrange=yrange,
                                   is.sen=sen,width=10.0,
                                   warrow="F",ave.x=ave.x,
                                   point.size=0.1)

    box(lwd=2.0)
          
    if ( i == 2 ) {
        axis(1,xaxp=xrange.axis,lwd=2.0,cex=1.0)
        mtext(outer=T,label.x,side=1,line=3.0,cex=1.0)
    }
        
    axis(2,yaxp=yrange.axis,lwd=2.0)
    mtext(outer=T,label.y,side=2,line=3.0,cex=1.0)
}

