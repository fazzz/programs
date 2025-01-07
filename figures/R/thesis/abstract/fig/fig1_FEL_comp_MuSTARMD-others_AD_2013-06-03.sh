#!~/bin/sh

proname=AD

dirOUT=~/thesis/abstract

parafile=~/calspa/refcalc/UmbSam/AD/para/para_s_vac_ff99SB_K=10_2012-11-12.sh
source ${parafile}

tiffwidth=800
tiffheight=600

nb=21

TLns=`expr ${TLbase} / 1000`
direqubase=${dirUMB}/s_UmbSam_vac_2012-11-12_${ff}/
dirpmf=${direqubase}${pname}/pmf
filenamepmf=${dirpmf}/pmf_UmbMD_vac_${TLns}ns.txt

paraRfile=${dirOUT}/graph/para_UmbSamp_AD_2Dpmf_2013-01-14.R

cat <<EOF > ${paraRfile}
fact.x <- 1
fact.y <- 1
fact.p <- 1.0/pi

dir    <- "${dirUMB}"
dirout <- "${dirOUT}"

title <- NULL

label.x <- expression(paste(""))
label.y <- expression(paste(""))

name.out <- paste(dirout,"/eps/","fig1_FEL_comp_MuSTARMD-others_AD_2013-06-03",sep='')
   
level <- seq(0.0,20.0,2.0)

file.name <- paste(name.out,'.eps',sep='')
postscript(file.name,width=3.2,height=5.8,horizontal=FALSE,onefile=FALSE,paper="special")

par(oma = c(0,0,0,0))

source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange_worwoaxis.R")
source("~/papers/CG-FG_TACCM_REMD/FillConBar.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap_wrange_worwoaxis.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConBar_wmai.R")

nf <- layout(matrix(c(1,4,2,4,3,4),3,2,byrow=TRUE),c(3,1),c(20,15,23))

par(cex.axis=1.44)
par(cex.lab=1.44)

EOF

#parafile_MuSTARMD=~/calspa/TACCM_CGAAREMD/AD/para/para_e_CG-FG_NH_2012-07-27_wrefd0.1_TZ=750_fq=10ps_99SB.sh 
#parafile_MuSTARMD=~/calspa/TACCM_CGAAREMD/AD/para/2012-11-02/para_e_CG-FG_NR=4_TZ=700_fq=10ps_99SB_KZmax=5000_weljd0.001_mZ=100.sh
parafile_MuSTARMD=~/calspa/TACCM_CGAAREMD/AD/para/2013-05-04/para_e_CG-FG_NR=8_2_TZ=750_fq=10ps_99SB_KZmax=5000_weljd0.001_mZ=100.sh

source ${parafile_MuSTARMD}
cat << eof >> ${paraRfile}
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
name <- paste(dirbase,"/tau=",tau_MuSTARMD,"/mZ=",mZ_MuSTARMD,"/TZ=",TZ_MuSTARMD,"/",pname_MuSTARMD,"/nEX=",numEX_MuSTARMD,"/fq=",fq_MuSTARMD,"/pmf/pmf_TAA=",TAA_MuSTARMD,"_TCG_",TCG_MuSTARMD,"_TZ_",TZ_MuSTARMD,"_",width_MuSTARMD,"_",KZAAo,"_",KZCGo,"_",AACG,"_4_2012-08-21",sep="") 
#name <- paste(dirbase,"/tau=",tau_MuSTARMD,"/mZ=",mZ_MuSTARMD,"/TZ=",TZ_MuSTARMD,"/",pname_MuSTARMD,"/nEX=",numEX_MuSTARMD,"/fq=",fq_MuSTARMD,"/pmf/pmf_pymbar_TAA=",TAA_MuSTARMD,"_TCG=",TCG_MuSTARMD,"_TZ=",TZ_MuSTARMD,"_KZAAo=",KZAAo,"_KZCGo=",KZCGo,"_",AACG,sep="") 

xrange <- c(-3.0,3.0,6)
xrange.axis <- c(-3.0,3.0,6)
yrange <- c(-3.0,3.0,6)
yrange.axis <- c(-2.0,3.0,5)

cat(name)

felwFillConMapwrangeworwoaxis(name,label.x=label.x,label.y=label.y,level=level,title=title,xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,mai=c(0.0,0.8,0.5,0.1),xaflag="f")

text(2.3,2.4,"(a)",cex=1.6)

eof

parafile_TREMD=~/calspa/refcalc/REMD/AD/para/para_s_REMD_vac_ff99SB_2012-07-24.sh
source ${parafile_TREMD}
cat <<EOF >> ${paraRfile}
pname_TREMD <-"${pname}"

numEX_TREMD <-"${numEX}"

TLbase_TREMD <-"${TLbase}"

ff_TREMD <-"${ff}"

dir0 <- "~/calspa/refcalc/REMD/AD"
dirbase <- paste(dir0,"/s_REVAC_2012-07-23_",ff_TREMD,sep="")
name <- paste(dirbase,"/",pname_TREMD,"/nEX=",numEX_TREMD,"/freq=",TLbase_TREMD,"ps","/pmf/pmf_ADv",sep="")

xrange <- c(-3.0,3.0,5)
xrange.axis <- c(-3.0,3.0,5)
yrange <- c(-3.0,3.0,5)
yrange.axis <- c(-2.0,3.0,5)

#label.x <- expression(paste(phi,"(radian)"))
label.y <- expression(paste(psi,"(radian)"))

felwFillConMapwrangeworwoaxis(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,mai=c(0.0,0.8,0.0,0.1),xaflag="f")

text(2.3,2.4,"(b)",cex=1.6)

EOF

parafile_TAMD=~/calspa/TACCM/AD/para/para_e_ff99SB_NH_2012-07-31_TZ=750.sh
source ${parafile_TAMD}
i=1
cat <<EOF >> ${paraRfile}
T_TAMD<-"${T}"

TB_TAMD<-"${TB[$i]}"

tau_TAMD<-"${tau}"

KZ_TAMD<-"${KZ}"

mZ_TAMD<-"${mZ}"

width_TAMD <- "0.3"

ff_TAMD <- "${ff}"

dir0 <- "/home/yamamori/calspa/TACCM/AD"
dirbase <- paste(dir0,"/e_TACCM_NH_2012-07-31_",ff_TAMD,sep="")
name <- paste(dirbase,"/tau=",tau_TAMD,"/TB=",TB_TAMD,"/KZ=",KZ_TAMD,"/mZ=",mZ_TAMD,"/pmf/pmf_T=",T_TAMD,".Zhist",sep="")

xrange <- c(-3.0,3.0,6)
xrange.axis <- c(-3.0,3.0,6)
yrange <- c(-3.0,3.0,6)
yrange.axis <- c(-3.0,3.0,6)

label.x <- expression(paste(phi,"(radian)"))
label.y <- ""

felwFillConMapwrangeworwoaxis(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,mai=c(0.8,0.8,0.0,0.1))

text(2.3,2.4,"(c)",cex=1.6)

EOF

cat <<EOF >> ${paraRfile}

felwFillConBarwmai(name,label.x=label.x,label.y=label.y,level=level,mai=c(0.8,0.5,0.5,0.2))

EOF

Rscript ${paraRfile}; echo ${paraRfile}
