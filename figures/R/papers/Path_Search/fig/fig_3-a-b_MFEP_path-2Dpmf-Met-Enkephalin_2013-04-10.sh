#!~/bin/sh

K=12

proname=MetEnk

minx=0.0
maxx=0.3
miny=0.0
maxy=0.3

xi=0.06
yi=0.15

xf=0.17
yf=0.06

source ~/calspa/MFEP/MetEnk/para/para_search_MFEPpmf2MGaussian_2013-03-26_bo_CG-FG_NH_2012-06-04.sh

height=5
widtho=1

dir=~/calspa/MFEP/MetEnk
dirf=~/calspa/MFEP/MetEnk/pmf2MGaussian_2013-04-10_bo_e_CG-FG_CMD_NH_2012-06-04

filenamepmf=${dirf}/pmf_TAA=${TAA}_TCG1_${TCG1}_TCG2_${TCG2}_TZ_${TZs}_${width}_${KZAAo}_${KZCG1o}_${KZCG2o}_${AACG}_bo10000ps_2

filenameMGaussian=${dirf}/pmf2MGaussian_TAA=${TAA}_TCG1=${TCG1}_TCG2=${TCG2}_${KZAAo}_${KZCG1o}_${KZCG2o}_K=${K}_@${xi},${yi}-@${xf},${yf}.txt

filenamepath=${dirf}/path_TAA=${TAA}_TCG1=${TCG1}_TCG2=${TCG2}_${KZAAo}_${KZCG1o}_${KZCG2o}_K=${K}_@${xi},${yi}-@${xf},${yf}.txt

parafile=~/papers/Path_Search/graph/fig_3-a-b_MFEP_path-2Dpmf-Met-Enkephalin_2013-04-10.R
cat <<EOF > ${parafile}
source("~/calspa/GOLMAA/SH3/R/filled.contour.mod.R")
source("~/calspa/GOLMAA/SH3/R/fel3.R")
source("~/Rspa/plRamachad3.R")

dir <- "~/papers/Path_Search"

title <- NULL

label.x <- expression(paste(d1))
label.y <- expression(paste(d2))

name.path <- paste("${filenamepath}",sep='')
name.out <- paste(dir,"/eps/","fig_3-a-b_MFEP_path-2Dpmf-Met-Enkephalin_2013-04-10",sep='')
   
level <- seq(0.0,${height},${widtho})

file.name <- paste(name.out,'.eps',sep='')
postscript(file.name,width=5.6,height=2.4,horizontal=FALSE,onefile=FALSE,paper="special")

xrange <- c(${minx},${maxx},3)
yrange <- c(${miny},${maxy},3)

par(oma = c(0,0,0,0) )
source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange.R")
source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange_wy.R")

source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMapwrangewxwy_wpath.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMapwrangewxwoy_wpath.R")

source("~/papers/CG-FG_TACCM_REMD/FillConBar.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConBar.R")

nf <- layout(matrix(c(1,2,3),1,3,byrow=TRUE),c(3.7,3,1.5),c(1))

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
