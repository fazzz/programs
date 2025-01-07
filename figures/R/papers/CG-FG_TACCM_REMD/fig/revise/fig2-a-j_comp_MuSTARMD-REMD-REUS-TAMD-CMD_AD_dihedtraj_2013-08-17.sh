#!~/bin/sh

parafileMuSTARMD=~/calspa/TACCM_CGAAREMD/AD/para/2012-11-02/para_e_CG-FG_NR=4_TZ=700_fq=10ps_99SB_KZmax=5000_weljd0.001_mZ=100.sh
#parafileTREMD=~/calspa/refcalc/REMD/AD/para/para_s_REMD_vac_ff99SB_2012-07-24.sh
parafileTREMD=~/calspa/refcalc/REMD/AD/para/para_s_REMD_vac_ff99SB_2013-08-17_TZ=300-1000.sh
parafileREUS=~/calspa/refcalc/REUS/AD/para/para_s_REUSMD_vac_ff99SB_Nbin=4x4_K=1_10ns_freq=10ps_2013-08-12.sh
parafileTAMD=~/calspa/TACCM/AD/para/para_e_ff99SB_NH_2012-07-31_TZ=750.sh
parafileCMD=~/calspa/refcalc/CMD/AD/para/para_s_CMD_T300_ff99SB_2013-02-28.sh

height=700
width=500

pointsize=0.1

dirMuSTARMD=~/calspa/TACCM_CGAAREMD/AD
dirTREMD=~/calspa/refcalc/REMD/AD
dirREUS=~/calspa/refcalc/REUS/AD
dirTAMD=~/calspa/TACCM/AD
dirCMD=~/calspa/refcalc/CMD/AD

dirOUT=~/papers/CG-FG_TACCM_REMD/

paraRfil=${dirOUT}/graph/graph_fig2-a-j_comp_MuSTARMD-REMD-REUS-TAMD-CMD_AD_dihedtraj_2013-08-17.R

source ${parafileMuSTARMD}

TAA=${TAA}
TCG=${TCG}
TZ=${TZs[1]}
tau=${tau[1]}
mZ=${mZ[1]}
pname=${pname}
TLbase=${TLbase}
numEX=${numEX}

cat <<EOF > ${paraRfil}
setwd("~/Rspa")
source("~/Rspa/ini.R")
source("~/Rspa/plmGeneral_wsrange_wdrange_wdiffx.R")

proname<-"AD"

TAA<-"${TAA}"
TCG<-"${TCG}"
TZ<-"${TZ}"

tau<-"${tau}"
mZ<-"${mZ}"

pname <-"${pname}"

numEX<-"${numEX}"
TLbase<-"${TLbase}"

name.title <- paste("${dirOUT}/eps/fig2-a-l_comp_MuSTARMD-REMD-REUS-TAMD-CMD_AD_dihedtraj_2013-08-17",sep="")

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

dir <- "${dirMuSTARMD}/e_CG-FG_NH_2012-07-27/"

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
                                   point.size=${pointsize})

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

EOF

source ${parafileTREMD}

ff=${ff}

Tbase=${Tbase}
TLbase=${TLbase}
numRE=${numRE}
pname=${pname}
numEX=${numEX}

cat <<EOF >> ${paraRfil}
TLbase<-"${TLbase}"
numEX<-"${numEX}"
numRE<-"${numRE}"
pname<-"${pname}"

filex.names<-NULL
filey.names<-NULL

yrange <- c(-1.0,1.0,4)
yrange.axis <- c(-1.0,0.5,3)

dir <- "${dirTREMD}/s_REVAC_2012-07-23_${ff}/"

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
                                   point.size=${pointsize})

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

EOF

source ${parafileREUS}

ff=${ff}

Tbase=${Tbase}
TLbase=${TLbase}
numREUS=${numREUS}
pname=${pname}
numEX=${numEX}

cat <<EOF >> ${paraRfil}
TLbase<-"${TLbase}"
numEX<-"${numEX}"
numREUS<-"${numREUS}"
pname<-"${pname}"

filex.names<-NULL
filey.names<-NULL

yrange <- c(-1.0,1.0,4)
yrange.axis <- c(-1.0,0.5,3)

dir <- "${dirREUS}/s_REUSVAC_2013-08-12_${ff}/"

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
                                   point.size=${pointsize})

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

EOF


source ${parafileTAMD}

ff=${ff}

tau=${tau}
TB=${TB[1]}
KZ=${KZ}
mZ=${mZ}
T=${T}

cat <<EOF >> ${paraRfil}
tau <-"${tau}"
TB <-"${TB}"
KZ <-"${KZ}"
mZ <-"${mZ}"
T <-"${T}"

filex.names<-NULL
filey.names<-NULL

yrange <- c(-1.0,1.0,4)
yrange.axis <- c(-1.0,0.5,3)

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

EOF

source ${parafileCMD}

ff=${ff}
T=${T}
TLbase=${TLbase}

cat <<EOF >> ${paraRfil}

T<-"${T}"
TLbase<-"${TLbase}"
numEX<-"${numEX}"
numRE<-"${numRE}"
pname<-"${pname}"

filex.names<-NULL
filey.names<-NULL

yrange <- c(-1.0,1.0,4)
yrange.axis <- c(-1.0,0.5,3)

dir <- "${dirCMD}/s_CVAC_2012-08-08_${ff}/"

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
                                   point.size=${pointsize})

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

EOF

Rscript ${paraRfil}
