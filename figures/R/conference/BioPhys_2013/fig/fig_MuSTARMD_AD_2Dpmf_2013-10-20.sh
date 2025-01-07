#!~/bin/sh

proname=AD
dirOUT=~/gakkai/BioPhys_2013

parafile=~/calspa/TACCM_CGAAREMD/AD/para/2013-08-28/para_e_CG-FG_NR=8_TZ=1000_fq=10ps_99SB_KZmax=1000_weljd0.001_mZ=100.sh

source ${parafile}

paraRfile=${dirOUT}/graph/para_MuSTARMD_AD_2Dpmf_2013-10-20.R

cat <<EOF > ${paraRfile}
TAA_MuSTARMD<-"${TAA}"
TCG_MuSTARMD<-"${TCG}"

TZ_MuSTARMD <- "${TZs[1]}"

tau_MuSTARMD <- "${tau[1]}"

mZ_MuSTARMD <-"${mZ[1]}"  

pname_MuSTARMD <- "${pname}"

numEX_MuSTARMD<-"${numEX}"

fq_MuSTARMD<-"${TLbase}"

TLbase_MuSTARMD<-"${TLbase}"

width_MuSTARMD <- "0.3"

dir0 <- "/home/yamamori/calspa/TACCM_CGAAREMD/AD"
dirbase <- paste(dir0,"/e_CG-FG_NH_2012-07-27",sep="")
AACG <- "CG"
KZAAo <- "0"
#KZCGo <- "1000"
KZCGo <- "5000"

name <- paste(dirbase,"/tau=",tau_MuSTARMD,"/mZ=",mZ_MuSTARMD,"/TZ=",TZ_MuSTARMD,"/",pname_MuSTARMD,"/nEX=",numEX_MuSTARMD,"/fq=",fq_MuSTARMD,"/pmf/pmf_pymbar_TAA=300_TCG=300_TZ=1000_KZAAo=1000_KZCGo=0_AA_2013-05-26_20_2",sep="") 

fact.x <- 180/pi
fact.y <- 180/pi
#fact.p <- 1.0/pi

dirout <- "${dirOUT}"

title <- NULL

label.x <- expression(paste(""))
label.y <- expression(paste(""))

name.out <- paste(dirout,"/eps/","fig_MuSTARMD_AD_2Dpmf_2013-10-20",sep='')
   
level <- seq(0.0,30.0,2.0)

xrange <- c(-180,180,6)
xrange.axis <- c(-180,180,6)
yrange <- c(-180,180,6)
yrange.axis <- c(-180,180,6)

file.name <- paste(name.out,'.eps',sep='')
postscript(file.name,width=4.0,height=2.5512,horizontal=FALSE,onefile=FALSE,paper="special")

par(oma = c(0,0,0,0) )
source("~/defense/fig/FillConMap.R")
source("~/defense/fig/FillConBar.R")
source("~/defense/fig/fel_FillConMap.R")
source("~/defense/fig/fel_FillConBar.R")

nf <- layout(matrix(c(1,2),1,2,byrow=TRUE),c(6.47,3.07),c(1))

par(cex.axis=1.2)
par(cex.lab=1.2)

cat(name)

label.x <- ""
label.y <- ""

felwFillConMap(name,label.x=label.x,label.y=label.y,level=level,
              title=title,kt="F/kBT",
              xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,
              mai=c(0.5,0.5,0.1,0.0))

felwFillConBar(name,label.x=label.x,label.y=label.y,level=level,
                   mai=c(1.348,0.75,0.1,0.75))

EOF

Rscript ${paraRfile}; echo ${paraRfile}
