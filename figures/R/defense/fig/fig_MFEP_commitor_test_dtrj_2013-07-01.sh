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

filenamepmf=~/calspa/TACCM_CGAAREMD/AD/e_CG-FG_NH_2012-07-27/tau=${tau[1]}/mZ=${mZ[1]}/TZ=${TZs[1]}/${pname}/nEX=${numEX}/fq=${TLbase}/pmf/pmf_TAA=${TAA}_TCG_${TCG}_TZ_${TZs[1]}_0.3_${KZAAo[$mode]}_${KZCGo[$mode]}_${MODES[$mode]}_4_2012-08-21

filenameMGaussian=${dir}/pmf2MGaussian_2012-10-10_bo_e_CG-FG_NH_2013-03-21/tau=${tau[1]}/mZ=${mZ[1]}/TZ=${TZs[1]}/${pname}/nEX=${numEX}/fq=${TLbase}/pmf2MGaussian_TAA=${TAA}_TCG=${TCG}_${KZAAo[$mode]}_${KZCGo[$mode]}_${MODES[$mode]}_K=${K}_wx${minx}-${maxx}_wy${miny}-${maxy}_dx=0.01_dy=0.01_2.txt

filenamepathbase=${dir}/committor_test_2013-03-24_bo_e_CG-FG_NH_2012-07-27/tau_${tau[1]}/mZ_${mZ[1]}/TZ_${TZs[1]}/${pname}/nEX_${numEX}/fq_${TLbase}/${crdTSL}/anl/${proname}v_T=${T}.dtrj

parafile=~/defense/graph/fig_MFEP_commitor_test_dtrj_2013-07-01.R

echo -n "name.path <- c(" > ${parafile}
for i in `seq 1 100`; do
    echo -n  "\"${filenamepathbase}_${i}\", " >> ${parafile}
done
echo  "sep='')" >> ${parafile}

cat <<EOF >> ${parafile}
source("~/calspa/GOLMAA/SH3/R/filled.contour.mod.R")
source("~/calspa/GOLMAA/SH3/R/fel3.R")
source("~/Rspa/plRamachad3.R")

dir <- "~/defense/"

title <- NULL

label.x <- expression(paste(phi))
label.y <- expression(paste(psi))

name.out <- paste(dir,"/eps/","fig_MFEP_commitor_test_dtrj_2013-07-01",sep='')
   
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

source("~/defense/fig/fel_wFillConMapwrangewxwy_wmdtrj.R")
source("~/defense/fig/fel_FillConBar.R")

nf <- layout(matrix(c(1,2),1,2,byrow=TRUE),c(6.47,3.07),c(1))

iro.base <- c(1,2,3,4,5,6,7,8,9,10)

iro <- rep(iro.base,10)

#par(mar =c(5,5,2,2) )
par(cex.axis=1.0)
par(cex.lab=1.5)

fact.x <- 180/pi
fact.y <- 180/pi
fact.p <- 180/pi

#name.pmf <- paste("${filenameMGaussian}",sep='')
name.pmf <- paste("${filenamepmf}",sep='')
cat(name.pmf," ",name.path,'\n')
name.path2 <- paste("${filenamepath2}",sep='')
fel_wFillConMapwrangewxwy_wmdtrj(name.pmf,name.path,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",mai=c(0.5,0.5,0.1,0.0),xrange=xrange,yrange=yrange,plot.axis="yes")

felwFillConBar(name.pmf,label.x=label.x,label.y=label.y,level=level,
                   mai=c(1.348,0.75,0.1,0.75))

EOF

Rscript ${parafile}


