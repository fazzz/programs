setwd("~/Rspa")
source("~/Rspa/ini.R")
source("~/Rspa/plmGeneral_wsrange_wdrange_wdiffx.R")

proname<-"AD"

T<-"300"
TLbase<-"80000"
numEX<-"100000"
numRE<-"8"
pname<-"f300t400"

name.title <- paste("/home/yamamori/papers/CG-FG_TACCM_REMD//tiff/fig5-c_CMD_AD_dihedtraj_2013-01-14",sep="")

leg.pos="topright"
is.leg <- 0

xrange <- c(0,10000,4)
xrange.axis <- c(2000,10000,2)

file.name <- paste(name.title,".tiff",sep="")
tiff(file.name,width=600,height=300)

par(mfrow=c(1,2))
par(mar=c(0.0,0.0,0.0,0.0))
par(oma=c(7.5,6.0,2.0,2.0))

filex.names<-NULL
filey.names<-NULL

yrange <- c(-3.15,3.15,10)
yrange.axis <- c(-3.0,3.0,6)

label.x <- expression(paste("MD Step (ps)"))
label.y <- "Dihedral Angle (radian)"

dir <- "/home/yamamori/calspa/refcalc/CMD/AD/s_CVAC_2012-08-08_ff99SB_2/"

hutosa <- NULL

for (i in 1:2) {
    sen<-c(0)
    senshu <- 1
    iro<-c(2)
    id.ys<-i

    dirbase.name <- paste(dir,"/",sep="")
        
    filex.names<-paste(dirbase.name,"anl/",proname,"v_T=",T,".dtrj_n",sep='')
    filey.names<-paste(dirbase.name,"anl/",proname,"v_T=",T,".dtrj",sep='')

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
                                   point.size=0.5)

    box(lwd=2.0)
          
    axis(1,xaxp=xrange.axis,lwd=2.0,cex=1.0)
    mtext(outer=T,label.x,side=1,line=3.0,cex=1.0)
        
    if ( i == 1 ) {
        axis(2,yaxp=yrange.axis,lwd=2.0,cex=1.0)
        mtext(label.y,side=2,line=3.0,cex=1.0)
    }
}

dev.off()

