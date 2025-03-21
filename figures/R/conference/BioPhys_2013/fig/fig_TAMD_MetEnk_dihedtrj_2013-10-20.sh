#!~/bin/sh

pointsize=0.1

dirOUT=~/gakkai/BioPhys_2013

paraRfil=${dirOUT}/graph/graph_TAMD_MetEnk_dihedtraj_2013-10-20.R

cat <<EOF > ${paraRfil}
setwd("~/Rspa")
source("~/Rspa/ini.R")
source("~/Rspa/plmGeneral_wsrange_wdrange_wdiffx.R")

proname<-"MetEnk"

name.title <- paste("${dirOUT}/eps/fig_TAMD_MetEnk_dihedtraj_2013-10-20",sep="")

leg.pos="topright"
is.leg <- 0

fact.y=1

file.name <- paste(name.title,".eps",sep="")
postscript(file.name,width=1.528,height=1.89728,horizontal=FALSE,onefile=FALSE,paper="special")

par(mfrow=c(2,1))
par(mar=c(0.0,0.0,0.0,0.0))
par(oma=c(0.5,0.5,0.1,0.1))

par(cex=1.0)

filex.names<-NULL
filey.names<-NULL

xrange <- c(0,10000,5)
xrange.axis <- c(0,10000,5)

yrange <- c(0.0,0.5,5)
yrange.axis <- c(0.0,0.5,5)

label.x <- expression(paste("time (ps)"))
label.y <- expression(paste("Dihedral Angle"))

hutosa <- NULL

for (i in 1:2) {
    sen<-c(0)
    senshu <- 1
    iro<-c(2)
    id.ys<-i+1

    filex.names <- "~/calspa/TACCM/MetEnk/e_TACCM_NH_2013-08-28_99SB_80ns/tau=1.0/TB=600/KZ=3000/mZ=100.00/anl/MetEnk_e_T=300.ddist_n"

    filey.names <- "~/calspa/TACCM/MetEnk/e_TACCM_NH_2013-08-28_99SB_80ns/tau=1.0/TB=600/KZ=3000/mZ=100.00/anl/MetEnk_e_T=300.ddist_n"

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
