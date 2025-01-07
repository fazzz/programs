#!~/bin/sh

proname=AD
dirUMB=~/calspa/refcalc/UmbSam/AD
dirOUT=~/papers/CG-FG_TACCM_REMD

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

#label.x <- expression(paste(phi,"(radian/",pi,")"))
#label.y <- expression(paste(psi,"(radian/",pi,")"))

label.x <- expression(paste(""))
label.y <- expression(paste(""))

name.out <- paste(dirout,"/eps/","fig3-a-g_comp_MuSTARMD-REMD-REUS-TAMD-CMD_AD_2Dpmf_2013-08-17",sep='')
   
#level <- seq(0.0,10.0,1.0)
level <- seq(0.0,20.0,2.0)

pdatanamedummy <- paste(dirout,"/fig/points_dummy.txt",sep='')
pdataname <- paste(dirout,"/fig/points.txt",sep='')
ldataname <- paste(dirout,"/fig/lines.txt",sep='')
ldataname2 <- paste(dirout,"/fig/lines2.txt",sep='')
ldataname3 <- paste(dirout,"/fig/lines3.txt",sep='')

xrange <- c(-1.0,1.0,4)
xrange.axis <- c(-1.0,1.0,4)
yrange <- c(-1.0,1.0,4)
yrange.axis <- c(-0.5,1.0,3)

file.name <- paste(name.out,'.eps',sep='')
postscript(file.name,width=5.3,height=6.0,horizontal=FALSE,onefile=FALSE,paper="special")

par(oma = c(0,0,0,0) )
source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange_worwoaxis.R")
source("~/papers/CG-FG_TACCM_REMD/FillConBar.R")

source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap_wrange_worwoaxis_wpoints.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConBar_wmai.R")

source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMapbox_wpoints_2.R")

nf <- layout(matrix(c(1,2,7,3,4,7,5,6,7),3,3,byrow=TRUE),c(22.5,16.5,8),c(20,15,23))

#par(cex.axis=1.0)
#par(cex.lab=1.0)

par(cex.axis=1.44)
par(cex.lab=1.44)

name <- paste("${filenamepmf}_2",sep='')

cat(name)

#label.x <- expression(paste(phi,"(radian/",pi,")"))
label.y <- ""

felwFillConMapwrangeworwoaxiswpoints(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,mai=c(0.0,0.8,0.5,0.0),xaflag="f",norm="T",pdata=pdatanamedummy,iro=c("white","white","white"))

text(0.75,0.8,"(a)",cex=1.6)

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

xrange <- c(-1.0,1.0,4)
xrange.axis <- c(-1.0,1.0,4)
yrange <- c(-1.0,1.0,4)
yrange.axis <- c(-1.0,0.5,3)

cat(name)

#label.x <- expression(paste(phi,"(radian/",pi,")"))
label.y <- ""

felwFillConMapwrangeworwoaxiswpoints(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,mai=c(0.0,0.0,0.5,0.1),xaflag="f",yaflag="f",norm="T",pdata=pdatanamedummy,iro=c("white","white","white"))

text(0.75,0.8,"(b)",cex=1.6)

eof

#parafile_TREMD=~/calspa/refcalc/REMD/AD/para/para_s_REMD_vac_ff99SB_2012-07-24.sh
parafile_TREMD=~/calspa/refcalc/REMD/AD/para/para_s_REMD_vac_ff99SB_2013-08-17_TZ=300-1000.sh
source ${parafile_TREMD}
cat <<EOF >> ${paraRfile}
pname_TREMD <-"${pname}"

numEX_TREMD <-"${numEX}"

TLbase_TREMD <-"${TLbase}"

ff_TREMD <-"${ff}"

dir0 <- "~/calspa/refcalc/REMD/AD"
dirbase <- paste(dir0,"/s_REVAC_2012-07-23_",ff_TREMD,sep="")
name <- paste(dirbase,"/",pname_TREMD,"/nEX=",numEX_TREMD,"/freq=",TLbase_TREMD,"ps","/pmf/pmf_py3",sep="")

xrange <- c(-1.0,1.0,4)
xrange.axis <- c(-1.0,1.0,4)
yrange <- c(-1.0,1.0,4)
yrange.axis <- c(-1.0,1.0,4)

#label.x <- expression(paste(phi,"(radian/",pi,")"))
#label.y <- expression(paste(psi,"(radian/",pi,")"))

cat(name)

felwFillConMapwrangeworwoaxiswpoints(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,mai=c(0.0,0.8,0.0,0.0),xaflag="f",norm="T",pdata=pdatanamedummy,iro=c("white","white","white"))

text(0.75,0.8,"(c)",cex=1.6)

EOF

parafile_REUS=~/calspa/refcalc/REUS/AD/para/para_s_REUSMD_vac_ff99SB_Nbin=4x4_K=1_10ns_freq=10ps_2013-08-12.sh
source ${parafile_REUS}
cat <<EOF >> ${paraRfile}
pname_REUS <-"${pname}"

numEX_REUS <-"${numEX}"

TLbase_REUS <-"${TLbase}"

ff_REUS <-"${ff}"

dir0 <- "~/calspa/refcalc/REUS/AD"
dirbase <- paste(dir0,"/s_REUSVAC_2013-08-12_",ff_REUS,sep="")
name <- paste(dirbase,"/",pname_REUS,"/nEX_",numEX_REUS,"/freq_",TLbase_REUS,"ps","/pmf/pmf_UmbMD_vac_nbx=40_nby=40.txt_2",sep="")

xrange <- c(-1.0,1.0,4)
xrange.axis <- c(-1.0,1.0,4)
yrange <- c(-1.0,1.0,4)
yrange.axis <- c(-1.0,1.0,4)


cat(name)

felwFillConMapwrangeworwoaxiswpoints(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,mai=c(0.0,0.0,0.0,0.1),xaflag="f",yaflag="f",norm="T",pdata=pdatanamedummy,iro=c("white","white","white"))

text(0.75,0.8,"(d)",cex=1.6)

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

xrange <- c(-1.0,0.5,3)
xrange.axis <- c(-1.0,0.5,3)
yrange <- c(-1.0,1.0,4)
yrange.axis <- c(-1.0,1.0,4)

#label.x <- expression(paste(phi,"(radian/",pi,")"))
#label.y <- ""

cat(name)

felwFillConMapwrangeworwoaxiswpoints(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,mai=c(0.8,0.8,0.0,0.0),norm="T",pdata=pdatanamedummy,iro=c("white","white","black"))

text(0.75,0.8,"(e)",cex=1.6)

EOF

parafile_CMD=~/calspa/refcalc/CMD/AD/para/para_s_CMD_T300_ff99SB_2012-07-24.sh
source ${parafile_CMD}
cat <<eof >> ${paraRfile}
pname_CMD<-"${pname}"

T_CMD<-"${T}"

ff_CMD<-"${ff}"

width_CMD <- "0.2"

dir0 <- "~/calspa/refcalc/CMD/AD"
dirbase <- paste(dir0,"/s_CVAC_2012-08-08_",ff_CMD,sep="")
name <- paste(dirbase,"/anl/pmf_ADv_T",T_CMD,"_",width_CMD,sep="")

xrange <- c(-1.0,1.0,4)
xrange.axis <- c(-1.0,1.0,4)
yrange <- c(-1.0,1.0,4)
yrange.axis <- c(-1.0,0.5,3)

#label.x <- expression(paste("                                 ",phi,"(rad"))
#label.y <- ""

cat(name)

felwFillConMapwrangeworwoaxisboxwpoints(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,mai=c(0.8,0.0,0.0,0.1),norm="T",yaflag="f",pdata=pdataname,ldata=ldataname,ldata2=ldataname2,ldata3=ldataname3,iro=c("pink","pink","pink","red","red","red"))

text(-0.6,0.8,"C5",cex=1.2,col="white")
text(-0.5,0.20,"C7eq",cex=1.2,col="white")
text(0.6,-0.2,"C7ax",cex=1.2)
text(0.75,0.8,"(f)",cex=1.6)

eof

cat <<EOF >> ${paraRfile}

felwFillConBarwmai(name,label.x=label.x,label.y=label.y,level=level,mai=c(0.8,0.5,0.5,0.2))

EOF

Rscript ${paraRfile}; echo ${paraRfile}
