#!~/bin/sh

#K=( dummy 3 4 5 6 )
K=( dummy  16 17 18 19 )
nK=5
mode=2
tiffwidth=900
tiffheight=200
height=10
width=1

proname=AD

minx=-2.5
maxx=2.5
miny=-2.5
maxy=2.5

source ~/calspa/MFEP/AD/para/para_e_CG-FG_NR=4_wrefd0.1_TZ=750_fq=10ps_99SB.sh

dir=~/calspa/MFEP/AD

filename=( dummy )

filename[1]=~/calspa/TACCM_CGAAREMD/AD/e_CG-FG_NH_2012-07-27/tau=${tau[1]}/mZ=${mZ[1]}/TZ=${TZs[1]}/${pname}/nEX=${numEX}/fq=${TLbase}/pmf/pmf_TAA=${TAA}_TCG_${TCG}_TZ_${TZs[1]}_0.3_${KZAAo[$mode]}_${KZCGo[$mode]}_${MODES[$mode]}_4_2012-08-21
filename[2]=${dir}/pmf2MGaussian_2012-10-10_bo_e_CG-FG_NH_2013-03-21/tau=${tau[1]}/mZ=${mZ[1]}/TZ=${TZs[1]}/${pname}/nEX=${numEX}/fq=${TLbase}/pmf2MGaussian_TAA=${TAA}_TCG=${TCG}_${KZAAo[$mode]}_${KZCGo[$mode]}_${MODES[$mode]}_K=${K[1]}_wx${minx}-${maxx}_wy${miny}-${maxy}.txt
filename[3]=${dir}/pmf2MGaussian_2012-10-10_bo_e_CG-FG_NH_2013-03-21/tau=${tau[1]}/mZ=${mZ[1]}/TZ=${TZs[1]}/${pname}/nEX=${numEX}/fq=${TLbase}/pmf2MGaussian_TAA=${TAA}_TCG=${TCG}_${KZAAo[$mode]}_${KZCGo[$mode]}_${MODES[$mode]}_K=${K[2]}_wx${minx}-${maxx}_wy${miny}-${maxy}.txt
filename[4]=${dir}/pmf2MGaussian_2012-10-10_bo_e_CG-FG_NH_2013-03-21/tau=${tau[1]}/mZ=${mZ[1]}/TZ=${TZs[1]}/${pname}/nEX=${numEX}/fq=${TLbase}/pmf2MGaussian_TAA=${TAA}_TCG=${TCG}_${KZAAo[$mode]}_${KZCGo[$mode]}_${MODES[$mode]}_K=${K[3]}_wx${minx}-${maxx}_wy${miny}-${maxy}.txt
filename[5]=${dir}/pmf2MGaussian_2012-10-10_bo_e_CG-FG_NH_2013-03-21/tau=${tau[1]}/mZ=${mZ[1]}/TZ=${TZs[1]}/${pname}/nEX=${numEX}/fq=${TLbase}/pmf2MGaussian_TAA=${TAA}_TCG=${TCG}_${KZAAo[$mode]}_${KZCGo[$mode]}_${MODES[$mode]}_K=${K[4]}_wx${minx}-${maxx}_wy${miny}-${maxy}.txt

parafile=~/papers/Path_Search/graph/fig_2DFEL-MixedGaussians_2012-10-24.R
cat <<EOF > ${parafile}
level <- seq(0,10,1)

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

name.out <- paste(dir,"/eps/","fig_2-a-b_2DFEL-MixedGaussians_2013-03-21",sep='')
file.name <- paste(name.out,'.eps',sep='')
postscript(file.name,width=5.6,height=2.0,horizontal=FALSE,onefile=FALSE,paper="special")

xrange <- c(-3,3,6)
yrange <- c(-3,3,6)

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
label.x <- expression(paste(phi))
label.y <- expression(paste(psi))

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
