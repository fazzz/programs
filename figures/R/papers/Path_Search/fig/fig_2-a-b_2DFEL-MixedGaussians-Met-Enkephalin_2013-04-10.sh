#!~/bin/sh

K=( dummy 4 6 12 13 14 )
nK=5

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
width=1

dir=~/calspa/MFEP/MetEnk
dirf=~/calspa/MFEP/MetEnk/pmf2MGaussian_2013-04-10_bo_e_CG-FG_CMD_NH_2012-06-04

filename=( dummy )

filename[1]=${dirf}/pmf2MGaussian_TAA=${TAA}_TCG1=${TCG1}_TCG2=${TCG2}_${KZAAo}_${KZCG1o}_${KZCG2o}_K=${K[1]}_@${xi},${yi}-@${xf},${yf}.txt
filename[2]=${dirf}/pmf2MGaussian_TAA=${TAA}_TCG1=${TCG1}_TCG2=${TCG2}_${KZAAo}_${KZCG1o}_${KZCG2o}_K=${K[2]}_@${xi},${yi}-@${xf},${yf}.txt
filename[3]=${dirf}/pmf2MGaussian_TAA=${TAA}_TCG1=${TCG1}_TCG2=${TCG2}_${KZAAo}_${KZCG1o}_${KZCG2o}_K=${K[3]}_@${xi},${yi}-@${xf},${yf}.txt
filename[4]=${dirf}/pmf2MGaussian_TAA=${TAA}_TCG1=${TCG1}_TCG2=${TCG2}_${KZAAo}_${KZCG1o}_${KZCG2o}_K=${K[4]}_@${xi},${yi}-@${xf},${yf}.txt
filename[5]=${dirf}/pmf2MGaussian_TAA=${TAA}_TCG1=${TCG1}_TCG2=${TCG2}_${KZAAo}_${KZCG1o}_${KZCG2o}_K=${K[5]}_@${xi},${yi}-@${xf},${yf}.txt

parafile=~/papers/Path_Search/graph/fig_2-a-b_2DFEL-MixedGaussians-Met-Enkephalin_2013-04-10.R
cat <<EOF > ${parafile}
level <- seq(0,${height},${width})

name.title <- NULL

MyColor <- function(n,alpha=1)
{			
    if ((n <- as.integer(n[1L])) > 0) {
        j <- n%/%3		
        k <- n%/%3		
        i <- n - j - k		
        c(if (i > 0) hsv(h = seq.int(from = 40/60, to = 25/60,
            length.out = i), alpha = alpha), if (j > 0) hsv(h = seq.int(from = 23/60,
            to = 11/60, length.out = j), alpha = alpha), if (k > 
            0) hsv(h = seq.int(from = 8/60, to = 0/60, length.out = k-1),
            alpha = alpha, s = seq.int(from = 1, to = 0.9, length.out = k-1),
            v = 1),hsv(0,0,1))
    }			
    else character(0L)	
}			

dir <- "~/papers/Path_Search"

title <- NULL

name <- paste("${filename}",sep='')

level <- seq(0.0,${height},${width})

name.out <- paste(dir,"/eps/","fig_2-a-b_2DFEL-MixedGaussians-Met-Enkephalin_2013-04-10",sep='')
file.name <- paste(name.out,'.eps',sep='')
postscript(file.name,width=5.6,height=2.0,horizontal=FALSE,onefile=FALSE,paper="special")

xrange <- c(${minx},${maxx},3)
yrange <- c(${miny},${maxy},3)

par(oma = c(0,0,0,0) )
source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange.R")
source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange_wy.R")
source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange_wox_wy.R")
source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange_woxyaxis.R")

source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap_wrange.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap_wrange_wy.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap_wrange_woxwy.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap_wrange_woxyaxis.R")

source("~/papers/CG-FG_TACCM_REMD/FillConBar.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConBar.R")

#nf <- layout(matrix(c(1,2,3,4,5,6),1,6,byrow=TRUE),c(3,3,3,3,3,1),c(1))
nf <- layout(matrix(c(1,2,3,4,5),1,5,byrow=TRUE),c(4.5,3,3,3,1.5),c(1))

#par(mar =c(5,5,2,2) )
par(cex.axis=1.0)
par(cex.lab=2.0)
label.x <- expression(paste(d1))
label.y <- expression(paste(d2))

#cat("${filename[1]}",'\n')
#felwFillConMap("${filename[1]}",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

cat("${filename[2]}",'\n')
felwFillConMapwrangewy("${filename[2]}",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",mar=c(5,5,2,0),xrange=xrange,yrange=yrange,plot.axis="yes")

cat("${filename[3]}",'\n')
felwFillConMapwrange("${filename[3]}",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",mar=c(5,0,2,0),xrange=xrange,yrange=yrange,plot.axis="yes")

cat("${filename[4]}",'\n')
felwFillConMapwrange("${filename[4]}",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",mar=c(5,0,2,0),xrange=xrange,yrange=yrange,plot.axis="yes")

cat("${filename[5]}",'\n')

felwFillConMapwrange("${filename[5]}",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",mar=c(5,0,2,0),xrange=xrange,yrange=yrange,plot.axis="yes")

par(mar=c(5,2,2,1))

felwFillConBar("${filename[1]}",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

EOF

Rscript ${parafile}
