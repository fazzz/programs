#!~/bin/sh

#K=6
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

proname=AD

#source ~/calspa/MFEP/AD/para/para_e_CG-FG_NR=4_wrefd0.1_TZ=350_fq=10ps_99SB_KZmax=1000.sh
source ~/calspa/MFEP/AD/para/para_e_CG-FG_NR=4_wrefd0.1_TZ=750_fq=10ps_99SB.sh

dir=~/calspa/MFEP/AD/

filenamepmf=~/calspa/TACCM_CGAAREMD/AD/e_CG-FG_NH_2012-07-27/tau=${tau[1]}/mZ=${mZ[1]}/TZ=${TZs[1]}/${pname}/nEX=${numEX}/fq=${TLbase}/pmf/pmf_TAA=${TAA}_TCG_${TCG}_TZ_${TZs[1]}_0.3_${KZAAo[$mode]}_${KZCGo[$mode]}_${MODES[$mode]}_4_2012-08-21_wrange_-2.0-2.0_-2.0-2.0

filenameMGaussian=${dir}/pmf2MGaussian_2012-10-10_bo_e_CG-FG_NH_2012-07-27/tau=${tau[1]}/mZ=${mZ[1]}/TZ=${TZs[1]}/${pname}/nEX=${numEX}/fq=${TLbase}/pmf2MGaussian_TAA=${TAA}_TCG=${TCG}_${KZAAo[$mode]}_${KZCGo[$mode]}_${MODES[$mode]}_K=${K}.txt

filenamepath=${dir}/pmf2MGaussian_2012-10-10_bo_e_CG-FG_NH_2012-07-27/tau=${tau[1]}/mZ=${mZ[1]}/TZ=${TZs[1]}/${pname}/nEX=${numEX}/fq=${TLbase}/path_TAA=${TAA}_TCG=${TCG}_${KZAAo[$mode]}_${KZCGo[$mode]}_${MODES[$mode]}_K=${K}_@${xi},${yi}-@${yi},${yf}.txt

parafile=~/Report/2012-10/graph/fig_MFEP_path-2Dpmf_2012-10-17.R
cat <<EOF > ${parafile}
source("~/calspa/GOLMAA/SH3/R/filled.contour.mod.R")
source("~/calspa/GOLMAA/SH3/R/fel3.R")
source("~/Rspa/plRamachad3.R")

dir <- "~/papers/Path_Search"

title <- NULL

label.x <- expression(paste(phi))
label.y <- expression(paste(psi))

name.path <- paste("${filenamepath}",sep='')
name.out <- paste(dir,"/eps/","fig_3-a-b_MFEP_path-2Dpmf_2013-03-19",sep='')
   
level <- seq(0.0,${height},${width})

file.name <- paste(name.out,'.eps',sep='')
postscript(file.name,width=7.0,height=2.8,horizontal=FALSE,onefile=FALSE,paper="special")

xrange <- c(-2,2,4)
yrange <- c(-2,2,4)

par(oma = c(0,0,0,0) )
source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange.R")
source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange_wy.R")

source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMapwrangewxwy_wpath.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMapwrangewxwoy_wpath.R")

source("~/papers/CG-FG_TACCM_REMD/FillConBar.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConBar.R")

nf <- layout(matrix(c(1,2,3),1,3,byrow=TRUE),c(3.7,3,1),c(1))

#par(mar =c(5,5,2,2) )
par(cex.axis=1.0)
par(cex.lab=2.0)

name.pmf <- paste("${filenameMGaussian}",sep='')
cat(name.pmf," ",name.path,'\n')
felwFillConMapwrangewxwywpath(name.pmf,name.path,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",mar=c(5,5,1,0),xrange=xrange,yrange=yrange,plot.axis="yes")

name.pmf <- paste("${filenamepmf}",sep='')
cat(name.pmf," ",name.path,'\n')
felwFillConMapwrangewxwoywpath(name.pmf,name.path,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,mar=c(5,0,1,0),plot.axis="yes")

par(mar=c(5,4,1,1))

felwFillConBar(name.pmf,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

EOF

Rscript ${parafile}


