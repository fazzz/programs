#!~/bin/sh

height=650
width=600

pointsize=0.1

dirOUT=~/papers/CG-FG_TACCM_REMD/

paraRfil=${dirOUT}/graph/graph_comp_MuSTARMD-REMD-REUS-TAMD-CMD_MetEnk_dihedangdistraj_2013-09-02.R

cat <<EOF > ${paraRfil}
setwd("~/Rspa")
source("~/Rspa/ini.R")
source("~/Rspa/plmGeneral_wsrange_wdrange_wdiffx.R")

name.title <- paste("${dirOUT}/eps/fig6-a-b_comp_MuSTARMD-CMD_MetEnk_dihedangdistraj_2013-09-02",sep="")

leg.pos="topright"
is.leg <- 0

xrange <- c(0,10000,5)

file.name <- paste(name.title,".eps",sep="")
#tiff(file.name,width=${width},height=${height})
postscript(file.name,width=6.1,height=7.7,horizontal=FALSE,onefile=FALSE,paper="special")

par(mfrow=c(5,2))
par(mar=c(0.0,0.0,0.0,0.0))
par(oma=c(4.0,4.0,1.0,1.0))
par(cex=1.0)

filex.names<-NULL
filey.names<-NULL

yrange <- c(0,0.5,5)

label.x <- expression(paste("MD Step (ps)"))
label.y <- "Dihedral Angle Distance"

hutosa <- NULL

for (i in 1:2) {
    sen<-c(0)
    senshu <- 1
    iro<-c(2)
    id.ys<-i+1

#    filex.names <- "~/calspa/TACCM_CGAAREMD/MetEnk/e_CG-FG_CMD_NH_2013-08-31//1FG2CG/tau=1.0/mZ=1000.00/ep=0.01/cutoff=4.7/TZ=600/TCG=400/2CG_2013-09-02_KAA=1000_KCG=1000-500/freq=1/fb=1.0/fa=1.0/ft=0.2/fc=1.0/fn=1.0/anl/MetEnk_Z_T=600_2.ddist_n"

#    filey.names <- "~/calspa/TACCM_CGAAREMD/MetEnk/e_CG-FG_CMD_NH_2013-08-31//1FG2CG/tau=1.0/mZ=1000.00/ep=0.01/cutoff=4.7/TZ=600/TCG=400/2CG_2013-09-02_KAA=1000_KCG=1000-500/freq=1/fb=1.0/fa=1.0/ft=0.2/fc=1.0/fn=1.0/anl/MetEnk_Z_T=600_2.ddist_n"

    filex.names <- "~/calspa/TACCM_CGAAREMD/MetEnk/e_CG-FG_CMD_NH_2013-08-31/1FG2CG/tau=1.0/mZ=1000.00/ep=0.01/cutoff=4.7/TZ=1000/TCG=400/2CG_2013-09-05_KAA=1250_KCG=1000-500_2/freq=1/fb=1.0/fa=1.0/ft=0.2/fc=1.0/fn=1.0/anl/MetEnk_Z_T=1000_1.ddist_n"

    filey.names <- "~/calspa/TACCM_CGAAREMD/MetEnk/e_CG-FG_CMD_NH_2013-08-31/1FG2CG/tau=1.0/mZ=1000.00/ep=0.01/cutoff=4.7/TZ=1000/TCG=400/2CG_2013-09-05_KAA=1250_KCG=1000-500_2/freq=1/fb=1.0/fa=1.0/ft=0.2/fc=1.0/fn=1.0/anl/MetEnk_Z_T=1000_1.ddist_n"

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
        txt <- "     (a)"
    }
    else {
        txt <- "     (b)"
    }
    text(200,0.37,txt)
          
    if ( i == 1 ) {
        yrange.axis <- c(0.1,0.4,3)
        axis(2,yaxp=yrange.axis,lwd=2.0)
#        mtext(label.y,side=2,line=3.0)
    }
}


for (i in 1:2) {
    sen<-c(0)
    senshu <- 1
    iro<-c(2)
    id.ys<-i+1

    filex.names <- "~/calspa/refcalc/REMD/MetEnk/s_RE_v_2012-07-25_SBL1/f300t600_80ns/nEX10000/frq1/anl/MetEnkv.CV_1_n"

    filey.names <- "~/calspa/refcalc/REMD/MetEnk/s_RE_v_2012-07-25_SBL1/f300t600_80ns/nEX10000/frq1/anl/MetEnkv.CV_1_n"

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
        txt <- "     (c)"
    }
    else {
        txt <- "     (d)"
    }
    text(200,0.37,txt)
          
    if ( i == 1 ) {
        yrange.axis <- c(0.1,0.4,3)
        axis(2,yaxp=yrange.axis,lwd=2.0)
#        mtext(label.y,side=2,line=3.0)
    }
}

#########################################################################################################################################
# for (i in 1:2) {														        #
#     sen<-c(0)															        #
#     senshu <- 1														        #
#     iro<-c(2)															        #
#     id.ys<-i+1														        #
# 																        #
#     filex.names <- "~/calspa/refcalc/REUS/MetEnk/s_REUSVAC_2013-08-30_ff99SB/UmbAt_Nbin_8_K_0.3/nEX_500/freq_10ps/anl/MetEnkv.CV_1_n" #
# 																        #
#     filey.names <- "~/calspa/refcalc/REUS/MetEnk/s_REUSVAC_2013-08-30_ff99SB/UmbAt_Nbin_8_K_0.3/nEX_500/freq_10ps/anl/MetEnkv.CV_1_n" #
# 																        #
#     hutosa <- 1.0														        #
#     ave.x <- 1														        #
#           															        #
#     plmGeneralwsrangewdrangewdiffx(datax.names=filex.names,									        #
#                                    datay.names=filey.names,									        #
#                                    sd.names=file.names,									        #
#                                    id.ys=id.ys,										        #
#                                    ids.ys=ids.ys,										        #
#                                    label.size=0.5,axis.size=2.0,								        #
#                                    iro=iro,axis.ft="F",is.header="F",								        #
#                                    sdiro=iro,											        #
#                                    xrange=xrange,yrange=yrange,								        #
#                                    sdyrange=yrange,										        #
#                                    is.sen=sen,width=10.0,									        #
#                                    warrow="F",ave.x=ave.x,									        #
#                                    point.size=${pointsize})									        #
# 																        #
#     box(lwd=2.0)														        #
#     if ( i == 1 ) {														        #
#         txt <- "     (e)"													        #
#     }																        #
#     else {															        #
#         txt <- "     (f)"													        #
#     }																        #
#     text(200,0.37,txt)													        #
#           															        #
#     if ( i == 1 ) {														        #
#         yrange.axis <- c(0.1,0.4,3)												        #
#         axis(2,yaxp=yrange.axis,lwd=2.0)											        #
# #        mtext(label.y,side=2,line=3.0)											        #
#     }																        #
# }																        #
#########################################################################################################################################

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
    if ( i == 1 ) {
        txt <- "     (e)"
    }
    else {
        txt <- "     (f)"
    }
    text(200,0.37,txt)
          
    if ( i == 1 ) {
        yrange.axis <- c(0.1,0.4,3)
        axis(2,yaxp=yrange.axis,lwd=2.0)
#        mtext(label.y,side=2,line=3.0)
    }
}

filex.names <- "~/calspa/refcalc/CMD/MetEnk/s_CVAC_2013-03-13_LM1/anl/MetEnk_T=300_fLM1.ddist_n"
              
filey.names <- "~/calspa/refcalc/CMD/MetEnk/s_CVAC_2013-03-13_LM1/anl/MetEnk_T=300_fLM1.ddist_1"

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
        txt <- "     (g)"
    }
    else {
        txt <- "     (h)"
    }
    text(200,0.37,txt)
          
    if ( i == 1 ) {
        yrange.axis <- c(0.1,0.4,3)
        axis(2,yaxp=yrange.axis,lwd=2.0)
#        mtext(label.y,side=2,line=3.0)
    }
}

filex.names <- "~/calspa/refcalc/CMD/MetEnk/s_CVAC_2013-03-13_LM2/anl/MetEnk_T=300_fLM2.ddist_n"
              
filey.names <- "~/calspa/refcalc/CMD/MetEnk/s_CVAC_2013-03-13_LM2/anl/MetEnk_T=300_fLM2.ddist_1"

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
        txt <- "     (i)"
    }
    else {
        txt <- "     (j)"
    }
    text(200,0.37,txt)

    if ( i == 1 ) {
        xrange.axis <- c(0,8000,4)
        axis(1,xaxp=xrange.axis,lwd=2.0)
        mtext(outer=T,label.x,side=1,line=3.0)
    }
    else {
        xrange.axis <- c(0,10000,5)
        axis(1,xaxp=xrange.axis,lwd=2.0)
    }
          
    if ( i == 1 ) {
        yrange.axis <- c(0.1,0.4,3)
        axis(2,yaxp=yrange.axis,lwd=2.0)
        mtext(outer=T,label.y,side=2,line=3.0)
    }
}

dev.off()

EOF

Rscript ${paraRfil}
