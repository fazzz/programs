#!~/bin/sh

height=650
width=600

pointsize=0.5

dirOUT=~/papers/CG-FG_TACCM_REMD/

paraRfil=${dirOUT}/graph/graph_comp_MuSTARMD-CMD_MetEnk_dihedangdistraj_2013-02-28.R

cat <<EOF > ${paraRfil}
setwd("~/Rspa")
source("~/Rspa/ini.R")
source("~/Rspa/plmGeneral_wsrange_wdrange_wdiffx.R")

name.title <- paste("${dirOUT}/tiff/fig6-a-b_comp_MuSTARMD-CMD_MetEnk_dihedangdistraj_2013-02-28",sep="")

leg.pos="topright"
is.leg <- 0

xrange <- c(0,10000,4)
xrange.axis <- c(2000,10000,2)

file.name <- paste(name.title,".tiff",sep="")
tiff(file.name,width=${width},height=${height})

par(mfrow=c(3,2))
par(mar=c(0.0,0.0,0.0,0.0))
par(oma=c(7.5,6.0,2.0,2.0))

filex.names<-NULL
filey.names<-NULL

yrange <- c(0,0.3,5)
yrange.axis <- c(0.06,0.3,4)

label.x <- expression(paste("MD Step (ps)"))
label.y <- "Dihedral Angle Distance"

hutosa <- NULL

for (i in 1:2) {
    sen<-c(0)
    senshu <- 1
    iro<-c(2)
    id.ys<-i

    filex.names <- "~/calspa/TACCM_CGAAREMD/MetEnk/e_CG-FG_CMD_NH_2012-06-04/1FG2CG/tau=1.0/mZ=50000.00/ep=0.2/cutoff=4.7/TZ=1400/2CG21_KAA=1000/freq=1/anl/MetEnk_Z_T=1400_3.ddist_n"

    filey.names <- "~/calspa/TACCM_CGAAREMD/MetEnk/e_CG-FG_CMD_NH_2012-06-04/1FG2CG/tau=1.0/mZ=50000.00/ep=0.2/cutoff=4.7/TZ=1400/2CG21_KAA=1000/freq=1/anl/MetEnk_Z_T=1400_3.ddist"

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
          
    if ( i == 1 ) {
        axis(2,yaxp=yrange.axis,lwd=2.0,cex=1.0)
        mtext(label.y,side=2,line=3.0,cex=1.0)
    }
}

filex.names <- "~/calspa/TACCM_CGAAREMD/MetEnk//e_CG-FG_NH_2012-05-14/FG_Amber/tau=1.0/anl/MetEnk_AA_T=300_fLM1.ddist_n"

filey.names <- "~/calspa/TACCM_CGAAREMD/MetEnk//e_CG-FG_NH_2012-05-14/FG_Amber/tau=1.0/anl/MetEnk_AA_T=300_fLM1.ddist"

hutosa <- NULL

for (i in 1:2) {
    sen<-c(0)
    senshu <- 1
    iro<-c(2)
    id.ys<-i

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
          
    if ( i == 1 ) {
        axis(2,yaxp=yrange.axis,lwd=2.0,cex=1.0)
        mtext(label.y,side=2,line=3.0,cex=1.0)
    }
}

filex.names <- "~/calspa/TACCM_CGAAREMD/MetEnk//e_CG-FG_NH_2012-05-14/FG_Amber/tau=1.0/anl/MetEnk_AA_T=300_fLM2.ddist_n"

filey.names <- "~/calspa/TACCM_CGAAREMD/MetEnk//e_CG-FG_NH_2012-05-14/FG_Amber/tau=1.0/anl/MetEnk_AA_T=300_fLM2.ddist"

hutosa <- NULL

for (i in 1:2) {
    sen<-c(0)
    senshu <- 1
    iro<-c(2)
    id.ys<-i

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
          
    axis(1,xaxp=xrange.axis,lwd=2.0,cex=1.0)
    mtext(outer=T,label.x,side=1,line=3.0,cex=1.0)
        
    if ( i == 1 ) {
        axis(2,yaxp=yrange.axis,lwd=2.0,cex=1.0)
        mtext(label.y,side=2,line=3.0,cex=1.0)
    }
}

dev.off()

EOF

Rscript ${paraRfil}
