#!~/bin/sh

K=18
mode=2
tiffwidth=700
tiffheight=300
height=10
width=1
xi=-1.4
yi=1.1
xf=1.1
yf=-0.8

minx=-4.0
maxx=4.0
miny=-4.0
maxy=4.0

proname=AD

source ~/calspa/MFEP/AD/para/para_e_CG-FG_NR=4_wrefd0.1_TZ=750_fq=10ps_99SB.sh

crdTSL=ADv_TS_-0.019--1.010.crd_3_1

dir=~/calspa/MFEP/AD/
dirout=~/defense

filenamepath=${dir}/pmf2MGaussian_2012-10-10_bo_e_CG-FG_NH_2013-03-21/tau=${tau[1]}/mZ=${mZ[1]}/TZ=${TZs[1]}/${pname}/nEX=${numEX}/fq=${TLbase}/path_TAA=${TAA}_TCG=${TCG}_${KZAAo[$mode]}_${KZCGo[$mode]}_${MODES[$mode]}_K=${K[4]}_@${xi},${yi}-@${yi},${yf}_wx${minx}-${maxx}_wy${miny}-${maxy}.txt

filenamepmf=~/calspa/TACCM_CGAAREMD/AD/e_CG-FG_NH_2012-07-27/tau=${tau[1]}/mZ=${mZ[1]}/TZ=${TZs[1]}/${pname}/nEX=${numEX}/fq=${TLbase}/pmf/pmf_TAA=${TAA}_TCG_${TCG}_TZ_${TZs[1]}_0.3_${KZAAo[$mode]}_${KZCGo[$mode]}_${MODES[$mode]}_4_2012-08-21

parafile=${dirout}/graph/fig_MFEP_pathon2DFEL_2013-07-09.R
cat <<EOF >> ${parafile}

dir <- "~/defense/"

title <- NULL

label.x <- expression(paste(phi))
label.y <- expression(paste(psi))

name.out <- paste(dir,"/eps/","fig_MFEP_pathon2DFEL_2013-07-01",sep='')
   
level <- seq(0.0,${height},${width})

file.name <- paste(name.out,'.eps',sep='')
postscript(file.name,width=4.1496,height=2.6472,horizontal=FALSE,onefile=FALSE,paper="special")

xrange <- c(-180,180,6)
xrange.axis <- c(-180,180,6)
yrange <- c(-180,180,6)
yrange.axis <- c(-180,180,6)

par(oma = c(0,0,0,0) )
source("~/defense/fig/FillConMap.R")
source("~/defense/fig/FillConBar.R")

source("~/defense/fig/fel_FillConMap_wpath.R")
source("~/defense/fig/fel_FillConBar.R")

nf <- layout(matrix(c(1,2),1,2,byrow=TRUE),c(6.47,3.07),c(1))

par(cex.axis=1.0)
par(cex.lab=1.5)

fact.x <- 180/pi
fact.y <- 180/pi
fact.p <- 180/pi

felwFillConMapwpath("${filenamepmf}","${filenamepath}",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",mai=c(0.5,0.5,0.1,0.0),xrange=xrange,yrange=yrange,plot.axis="yes")

felwFillConBar("${filenamepmf}",label.x=label.x,label.y=label.y,level=level,mai=c(1.348,0.75,0.1,0.75))

EOF

Rscript ${parafile}


