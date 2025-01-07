#!~/bin/sh

parafileTAMD=~/calspa/TACCM/AD/para/para_e_ff99SB_NH_2012-07-31_TZ=750.sh

pointsize=0.1

dirTAMD=~/calspa/TACCM/AD
dirOUT=~/defense/

paraRfil=${dirOUT}/graph/graph_MuSTARMD_AD_dihedtraj_2013-07-01.R

source ${parafileTAMD}

ff=${ff}

tau=${tau}
TB=${TB[1]}
KZ=${KZ}
mZ=${mZ}
T=${T}

cat <<EOF > ${paraRfil}
setwd("~/Rspa")
source("~/Rspa/ini.R")
source("~/Rspa/plmGeneral_wsrange_wdrange_wdiffx.R")

proname<-"AD"

tau <-"${tau}"
TB <-"${TB}"
KZ <-"${KZ}"
mZ <-"${mZ}"
T <-"${T}"

name.title <- paste("${dirOUT}/eps/fig_TAMD_AD_dihedtraj_2013-07-01",sep="")

leg.pos="topright"
is.leg <- 0

fact.y=1/pi

file.name <- paste(name.title,".eps",sep="")
postscript(file.name,width=1.5843,height=2.3716,horizontal=FALSE,onefile=FALSE,paper="special")

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

dir <- "${dirTAMD}/e_TACCM_NH_2012-07-31_${ff}/"

hutosa <- NULL

for (i in 1:2) {
    sen<-c(0)
    senshu <- 1
    iro<-c(2)
    id.ys<-i

    dirbase.name <- paste(dir,"/tau=",tau,"/TB=",TB,"/KZ=",KZ,"/mZ=",mZ,sep="")
        
    filex.names<-paste(dirbase.name,"/",proname,"_e_T=",T,".Theta_thin_n",sep='')
    filey.names<-paste(dirbase.name,"/",proname,"_e_T=",T,".Theta_thin",sep='')

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
                                   point.size=${pointsize})

    box(lwd=2.0)
          
    if ( i == 2 ) {
        axis(1,xaxp=xrange.axis,lwd=2.0,cex=1.0)
        mtext(outer=T,label.x,side=1,line=3.0,cex=1.0)
    }
        
    axis(2,yaxp=yrange.axis,lwd=2.0)
    mtext(outer=T,label.y,side=2,line=3.0,cex=1.0)
}

EOF

Rscript ${paraRfil}
