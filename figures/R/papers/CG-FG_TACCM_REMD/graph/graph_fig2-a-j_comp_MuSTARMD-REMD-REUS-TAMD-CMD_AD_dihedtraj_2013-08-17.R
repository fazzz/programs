setwd("~/Rspa")
source("~/Rspa/ini.R")
source("~/Rspa/plmGeneral_wsrange_wdrange_wdiffx.R")

proname<-"AD"

TAA<-"300"
TCG<-"300"
TZ<-"700"

tau<-"1.0"
mZ<-"100.00"

pname <-"SB_KZMAX=5000_NR4_woeljd0.001"

numEX<-"1000"
TLbase<-"10"

name.title <- paste("/home/yamamori/papers/CG-FG_TACCM_REMD//eps/fig2-a-l_comp_MuSTARMD-REMD-REUS-TAMD-CMD_AD_dihedtraj_2013-08-17",sep="")

leg.pos="topright"
is.leg <- 0

fact.y=1/pi

file.name <- paste(name.title,".eps",sep="")
postscript(file.name,width=6.1,height=7.5,horizontal=FALSE,onefile=FALSE,paper="special")

par(mfrow=c(5,2))
par(mar=c(0.0,0.0,0.0,0.0))
par(oma=c(4.0,4.0,1.0,1.0))

par(cex=1.0)

filex.names<-NULL
filey.names<-NULL

xrange <- c(0,10000,5)

yrange <- c(-1.0,1.0,4)
yrange.axis <- c(-1.0,1.0,4)

#################################################################
# label.x <- expression(paste("Time (ps)"))		        #
# label.y <- expression(paste("Dihedral Angle (radian",pi,")")) #
#################################################################

label.x <- " "
label.y <- " "

dir <- "/home/yamamori/calspa/TACCM_CGAAREMD/AD/e_CG-FG_NH_2012-07-27/"

hutosa <- NULL

for (i in 1:2) {
    sen<-c(0)
    senshu <- 1
    iro<-c(2)
    id.ys<-i

    dirbase.name <- paste(dir,"/tau=",tau,"/mZ=",mZ,"/TZ=",TZ,"/",pname,"/nEX=",numEX,"/fq=",TLbase,"/",sep="")
        
#    filex.names<-paste(dirbase.name,proname,"_AA_T=",TAA,".out_1",sep='')
    filex.names<-paste(dirbase.name,proname,"_AA.Theta_1_n",sep='')
    filey.names<-paste(dirbase.name,proname,"_AA.Theta_1",sep='')

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
        text(200,0.85,txt)
    }
    else {
        txt <- "     (b)"
        text(200,0.85,txt,font=2,col="white")
#        text(200,0.85,txt,font=2,col="black")
    }
          
#    axis(1,xaxp=xrange.axis,lwd=2.0,cex=1.0)
#    mtext(outer=T,label.x,side=1,line=3.0,cex=1.0)
        
    if ( i == 1 ) {
        axis(2,yaxp=yrange.axis,lwd=2.0)
#        mtext(label.y,side=2,line=3.0,cex=1.0)
    }
}

TLbase<-"1"
numEX<-"100000"
numRE<-"8"
pname<-"f300t1000"

filex.names<-NULL
filey.names<-NULL

yrange <- c(-1.0,1.0,4)
yrange.axis <- c(-1.0,0.5,3)

dir <- "/home/yamamori/calspa/refcalc/REMD/AD/s_REVAC_2012-07-23_ff99SB/"

hutosa <- NULL

for (i in 1:2) {
    sen<-c(0)
    senshu <- 1
    iro<-c(2)
    id.ys<-i

    dirbase.name <- paste(dir,"/",pname,"/nEX=",numEX,"/freq=",TLbase,"ps/",sep="")
        
    filex.names<-paste(dirbase.name,"anl/",proname,"g.CV_1_n",sep='')
    filey.names<-paste(dirbase.name,"anl/",proname,"g.CV_1",sep='')

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
        text(200,0.85,txt)
    }
    else {
        txt <- "     (d)"
        text(200,0.85,txt,font=2,col="white")
#        text(200,0.85,txt,font=2,col="black")
    }
          
    if ( i == 1 ) {
        axis(2,yaxp=yrange.axis,lwd=2.0)
#        mtext(label.y,side=2,line=3.0,cex=1.0)
    }
}

TLbase<-"10"
numEX<-"500"
numREUS<-"16"
pname<-"UmbAt_Nbin_4x4_K_1"

filex.names<-NULL
filey.names<-NULL

yrange <- c(-1.0,1.0,4)
yrange.axis <- c(-1.0,0.5,3)

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
    if ( i == 1 ) {
        txt <- "     (e)"
        text(200,0.85,txt)
    }
    else {
        txt <- "     (f)"
        text(200,0.85,txt,font=2,col="black")
#        text(200,0.85,txt,font=2,col="black")
    }
          
    if ( i == 1 ) {
        axis(2,yaxp=yrange.axis,lwd=2.0)
#        mtext(label.y,side=2,line=3.0,cex=1.0)
    }
}

tau <-"1.0"
TB <-"750"
KZ <-"1000"
mZ <-"50000.00"
T <-"300"

filex.names<-NULL
filey.names<-NULL

yrange <- c(-1.0,1.0,4)
yrange.axis <- c(-1.0,0.5,3)

dir <- "/home/yamamori/calspa/TACCM/AD/e_TACCM_NH_2012-07-31_99SB/"

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
                                   point.size=0.1)

    box(lwd=2.0)
    if ( i == 1 ) {
        txt <- "     (g)"
        text(200,0.85,txt)
    }
    else {
        txt <- "     (h)"
        text(120,0.85,txt)
    }
          
    if ( i == 1 ) {
        axis(2,yaxp=yrange.axis,lwd=2.0)
#        mtext(label.y,side=2,line=3.0,cex=1.0)
    }
}


T<-"300"
TLbase<-"800000"
numEX<-"500"
numRE<-"8"
pname<-"UmbAt_Nbin_4x4_K_1"

filex.names<-NULL
filey.names<-NULL

yrange <- c(-1.0,1.0,4)
yrange.axis <- c(-1.0,0.5,3)

dir <- "/home/yamamori/calspa/refcalc/CMD/AD/s_CVAC_2012-08-08_ff99SB_3/"

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
                                   point.size=0.1)

    box(lwd=2.0)
    if ( i == 1 ) {
        txt <- "     (i)"
        text(200,0.85,txt)
    }
    else {
        txt <- "     (j)"
        text(200,0.85,txt)
    }

    if ( i == 1 ) {
        xrange.axis <- c(0,8000,4)
    }
    else {
        xrange.axis <- c(0,10000,5)
    }
          
    axis(1,xaxp=xrange.axis,lwd=2.0)
    mtext(outer=T,label.x,side=1,line=3.0,cex=1.0)
        
    if ( i == 1 ) {
        axis(2,yaxp=yrange.axis,lwd=2.0)
        mtext(outer=T,label.y,side=2,line=3.0,cex=1.0)
    }

}
#dev.off()

