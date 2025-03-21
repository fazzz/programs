#!~/bin/sh

parafileMuSTARMD=~/calspa/TACCM_CGAAREMD/AD/para/2012-11-02/para_e_CG-FG_NR=4_TZ=700_fq=10ps_99SB_KZmax=5000_weljd0.001_mZ=100.sh
parafileTREMD=~/calspa/refcalc/REMD/AD/para/para_s_REMD_vac_ff99SB_2012-07-24.sh
#parafileCMD=~/calspa/refcalc/CMD/AD/para/para_s_CMD_T300_ff99SB_2012-07-24.sh
parafileCMD=~/calspa/refcalc/CMD/AD/para/para_s_CMD_T300_ff99SB_2013-02-28.sh

height=700
width=500

pointsize=0.5

dirMuSTARMD=~/calspa/TACCM_CGAAREMD/AD
dirTREMD=~/calspa/refcalc/REMD/AD
dirCMD=~/calspa/refcalc/CMD/AD

dirOUT=~/papers/CG-FG_TACCM_REMD/

paraRfil=${dirOUT}/graph/graph_MuSTARMD-TREMD-CMD_AD_dihedtraj_2013-02-28.R

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

name.title <- paste("${dirOUT}/eps/fig2-a-b-c_comp_MuSTARMD-TREMD-CMD_AD_dihedtraj_2013-02-28",sep="")

leg.pos="topright"
is.leg <- 0

xrange <- c(0,10000,4)
xrange.axis <- c(2000,10000,2)

file.name <- paste(name.title,".eps",sep="")
#tiff(file.name,width=${width},height=${height})
postscript(file.name,width=3.5,height=5.2,horizontal=FALSE,onefile=FALSE,paper="special")

par(mfrow=c(4,2))
par(mar=c(0.0,0.0,0.0,0.0))
par(oma=c(7.5,6.0,2.0,2.0))

filex.names<-NULL
filey.names<-NULL

yrange <- c(-3.15,3.15,10)
yrange.axis <- c(-3.0,3.0,6)

label.x <- expression(paste("MD Step (ps)"))
label.y <- "Dihedral Angle (radian)"

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
          
#    axis(1,xaxp=xrange.axis,lwd=2.0,cex=1.0)
#    mtext(outer=T,label.x,side=1,line=3.0,cex=1.0)
        
    if ( i == 1 ) {
        axis(2,yaxp=yrange.axis,lwd=2.0,cex=0.5)
        mtext(label.y,side=2,line=3.0,cex=0.5)
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
        axis(2,yaxp=yrange.axis,lwd=2.0,cex=0.5)
        mtext(label.y,side=2,line=3.0,cex=0.5)
    }
}

TLbase<-"${TLbase}"
numEX<-"${numEX}"
numRE<-"${numRE}"
pname<-"${pname}"

filex.names<-NULL
filey.names<-NULL

dir <- "${dirTREMD}/s_REVAC_2012-07-23_${ff}/"

hutosa <- NULL

for (i in 1:2) {
    sen<-c(0)
    senshu <- 1
    iro<-c(2)
    id.ys<-i

    dirbase.name <- paste(dir,"/",pname,"/nEX=",numEX,"/freq=",TLbase,"ps/",sep="")
        
    filex.names<-paste(dirbase.name,"anl/",proname,"g.CV_replica6_n",sep='')
    filey.names<-paste(dirbase.name,"anl/",proname,"g.CV_replica6",sep='')

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
        axis(2,yaxp=yrange.axis,lwd=2.0,cex=0.5)
        mtext(label.y,side=2,line=3.0,cex=0.5)
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

yrange <- c(-3.15,3.15,10)
yrange.axis <- c(-3.0,3.0,6)

label.x <- expression(paste("MD Step (ps)"))
label.y <- "Dihedral Angle (radian)"

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
          
    axis(1,xaxp=xrange.axis,lwd=2.0,cex=0.4)
    mtext(outer=T,label.x,side=1,line=3.0,cex=0.4)
        
    if ( i == 1 ) {
        axis(2,yaxp=yrange.axis,lwd=2.0,cex=0.5)
        mtext(label.y,side=2,line=3.0,cex=0.5)
    }

}
#dev.off()

EOF

Rscript ${paraRfil}
