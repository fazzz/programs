setwd("~/Rspa")
source("~/Rspa/ini.R")
source("~/Rspa/plmGeneral_wsrange_wdrange_wdiffx.R")

name.title <- paste("/home/yamamori/papers/CG-FG_TACCM_REMD//eps/fig6-a-b_comp_MuSTARMD-CMD_MetEnk_dihedangdistraj_2013-06-04",sep="")

leg.pos="topright"
is.leg <- 0

xrange <- c(0,10000,5)

file.name <- paste(name.title,".eps",sep="")
#tiff(file.name,width=600,height=650)
#postscript(file.name,width=7.0,height=6.3,horizontal=FALSE,onefile=FALSE,paper="special")
postscript(file.name,width=6.1,height=6.3,horizontal=FALSE,onefile=FALSE,paper="special")

par(mfrow=c(3,2))
par(mar=c(0.0,0.0,0.0,0.0))
#par(oma=c(7.5,6.0,2.0,2.0))
#par(oma=c(5.0,5.0,2.0,2.0))
#par(oma=c(4.0,4.0,1.0,3.0))
par(oma=c(4.0,4.0,1.0,1.0))
par(cex=1.0)

filex.names<-NULL
filey.names<-NULL

yrange <- c(0,0.4,4)

label.x <- expression(paste("MD Step (ps)"))
label.y <- "Dihedral Angle Distance"

hutosa <- NULL

for (i in 1:2) {
    sen<-c(0)
    senshu <- 1
    iro<-c(2)
    id.ys<-i

    filex.names <- "~/calspa/TACCM_CGAAREMD/MetEnk/e_CG-FG_CMD_NH_2012-06-04/1FG2CG/tau=1.0/mZ=1000.00/ep=0.2/cutoff=4.7/TZ=2000/2CG_2013-05-31_KZ=2000/freq=1/anl/MetEnk_Z_T=2000_1.ddist_n"

    filey.names <- "~/calspa/TACCM_CGAAREMD/MetEnk/e_CG-FG_CMD_NH_2012-06-04/1FG2CG/tau=1.0/mZ=1000.00/ep=0.2/cutoff=4.7/TZ=2000/2CG_2013-05-31_KZ=2000/freq=1/anl/MetEnk_Z_T=2000_1.ddist"

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

#filex.names <- "~/calspa/TACCM_CGAAREMD/MetEnk//e_CG-FG_NH_2012-05-14/FG_Amber/tau=1.0/anl/MetEnk_AA_T=300_fLM1.ddist_n"

#filey.names <- "~/calspa/TACCM_CGAAREMD/MetEnk//e_CG-FG_NH_2012-05-14/FG_Amber/tau=1.0/anl/MetEnk_AA_T=300_fLM1.ddist"

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
                                   point.size=0.1)

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

#filex.names <- "~/calspa/TACCM_CGAAREMD/MetEnk//e_CG-FG_NH_2012-05-14/FG_Amber/tau=1.0/anl/MetEnk_AA_T=300_fLM2.ddist_n"

#filey.names <- "~/calspa/TACCM_CGAAREMD/MetEnk//e_CG-FG_NH_2012-05-14/FG_Amber/tau=1.0/anl/MetEnk_AA_T=300_fLM2.ddist"

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
                                   point.size=0.1)

    box(lwd=2.0)
    if ( i == 1 ) {
        txt <- "     (e)"
    }
    else {
        txt <- "     (f)"
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
        yrange.axis <- c(0.0,0.4,4)
        axis(2,yaxp=yrange.axis,lwd=2.0)
        mtext(outer=T,label.y,side=2,line=3.0)
    }
}

dev.off()

