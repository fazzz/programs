#!~/bin/sh

K=17
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

source ~/calspa/MFEP/AD/para/para_e_CG-FG_NR=8_2_TZ=750_fq=10ps_99SB_KZmax=5000_weljd0.001_mZ=100_2013-06-05.sh

dir=~/calspa/MFEP/AD/

filenamepmf=~/calspa/TACCM_CGAAREMD/AD/e_CG-FG_NH_2012-07-27/tau=${tau[1]}/mZ=${mZ[1]}/TZ=${TZs[1]}/${pname}/nEX=${numEX}/fq=${TLbase}/pmf/pmf_TAA=${TAA}_TCG_${TCG}_TZ_${TZs[1]}_0.3_${KZAAo[$mode]}_${KZCGo[$mode]}_${MODES[$mode]}_4_2012-08-21_wrange_${minx}-${maxx}_${miny}-${maxy}

filenameMGaussian=${dir}/pmf2MGaussian_2012-10-10_bo_e_CG-FG_NH_2013-05-04/tau=${tau[1]}/mZ=${mZ[1]}/TZ=${TZs[1]}/${pname}/nEX=${numEX}/fq=${TLbase}/pmf2MGaussian_TAA=${TAA}_TCG=${TCG}_${KZAAo[$mode]}_${KZCGo[$mode]}_${MODES[$mode]}_K=${K}_wx${minx}-${maxx}_wy${miny}-${maxy}_dx=0.01_dy=0.01.txt

filenamepath2=${dir}/pmf2MGaussian_2012-10-10_bo_e_CG-FG_NH_2013-05-04/tau=${tau[1]}/mZ=${mZ[1]}/TZ=${TZs[1]}/${pname}/nEX=${numEX}/fq=${TLbase}/path_TAA=${TAA}_TCG=${TCG}_${KZAAo[$mode]}_${KZCGo[$mode]}_${MODES[$mode]}_K=${K}_@${xi},${yi}-@${yi},${yf}_wx${minx}-${maxx}_wy${miny}-${maxy}.txt


for m in `seq 1 ${ncrdTS}`; do
    filenamepathbase=${dir}/committor_test_2013-03-24_bo_e_CG-FG_NH_2012-07-27/tau_${tau[1]}/mZ_${mZ[1]}/TZ_${TZs[1]}/${pname}/nEX_${numEX}/fq_${TLbase}/${crdTS[$m]}/anl/${proname}v_T=${T}.dtrj

    parafile=~/papers/Path_Search/graph/fig_5_committor_test_2013-03-24.R

    echo -n "name.path <- c(" > ${parafile}
    for i in `seq 1 100`; do
	echo -n  "\"${filenamepathbase}_${i}\", " >> ${parafile}
    done
    echo  "sep='')" >> ${parafile}

    cat <<EOF >> ${parafile}
source("~/calspa/GOLMAA/SH3/R/filled.contour.mod.R")
source("~/calspa/GOLMAA/SH3/R/fel3.R")
source("~/Rspa/plRamachad3.R")

dir <- "~/papers/Path_Search"

title <- NULL

label.x <- expression(paste(phi))
label.y <- expression(paste(psi))

name.out <- paste(dir,"/eps/","fig_5_committor_test_2013-06-07_${crdTS[$m]}",sep='')
   
level <- seq(0.0,${height},${width})

file.name <- paste(name.out,'.eps',sep='')
postscript(file.name,width=4.2,height=3.4,horizontal=FALSE,onefile=FALSE,paper="special")

xrange <- c(${minx},${maxx},5)
yrange <- c(${miny},${maxx},5)

par(oma = c(0,0,0,0) )
source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange.R")
source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange_wy.R")

source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMapwrangewxwy_wmpath2.R")

source("~/papers/CG-FG_TACCM_REMD/FillConBar.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConBar.R")

nf <- layout(matrix(c(1,2),1,2,byrow=TRUE),c(3.7,1.5),c(1))

iro.base <- c(1,2,3,4,5,6,7,8,9,10)

iro <- rep(iro.base,10)

#par(mar =c(5,5,2,2) )
par(cex.axis=1.0)
par(cex.lab=1.5)

name.pmf <- paste("${filenameMGaussian}",sep='')
cat(name.pmf," ",name.path,'\n')
name.path2 <- paste("${filenamepath2}",sep='')
felwFillConMapwrangewxwywmpath2(name.pmf,name.path,name.path2,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",mar=c(5,5,1,0),xrange=xrange,yrange=yrange,plot.axis="yes")

par(mar=c(5,4,1,1))

felwFillConBar(name.pmf,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

EOF

    Rscript ${parafile}

done
