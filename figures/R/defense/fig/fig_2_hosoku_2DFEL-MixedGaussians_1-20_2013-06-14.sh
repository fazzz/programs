#!~/bin/sh

mode=2
tiffwidth=900
tiffheight=200
height=10
width=1
#height=20
#width=2

proname=AD

xi=-1.4
yi=1.1
xf=1.1
yf=-0.8

minx=-4.0
maxx=4.0
miny=-4.0
maxy=4.0

source ~/calspa/MFEP/AD/para/para_e_CG-FG_NR=8_2_TZ=750_fq=10ps_99SB_KZmax=5000_weljd0.001_mZ=100_2013-06-05.sh

dir=~/calspa/MFEP/AD

filename=( dummy )

filename[1]=~/calspa/TACCM_CGAAREMD/AD/e_CG-FG_NH_2012-07-27/tau=${tau[1]}/mZ=${mZ[1]}/TZ=${TZs[1]}/${pname}/nEX=${numEX}/fq=${TLbase}/pmf/pmf_TAA=${TAA}_TCG_${TCG}_TZ_${TZs[1]}_0.3_${KZAAo[$mode]}_${KZCGo[$mode]}_${MODES[$mode]}_4_2012-08-21

for i in `seq 2 20`; do
    filename[$i]=${dir}/pmf2MGaussian_2012-10-10_bo_e_CG-FG_NH_2013-05-04/tau=${tau[1]}/mZ=${mZ[1]}/TZ=${TZs[1]}/${pname}/nEX=${numEX}/fq=${TLbase}/pmf2MGaussian_TAA=${TAA}_TCG=${TCG}_${KZAAo[$mode]}_${KZCGo[$mode]}_${MODES[$mode]}_K=${i}_wx${minx}-${maxx}_wy${miny}-${maxy}_dx=0.01_dy=0.01_4.txt
done

parafile=~/papers/Path_Search/graph/fig_2_hosoku_2DFEL-MixedGaussians_1-20_2013-06-14.R
cat <<EOF > ${parafile}
fact.x <- 1
fact.y <- 1
fact.p <- 1.0/pi

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

name.out <- paste(dir,"/eps/","fig_2_hosoku_2DFEL-MixedGaussians_1-20_2013-06-14",sep='')
file.name <- paste(name.out,'.eps',sep='')
#postscript(file.name,width=6.4,height=8.0,horizontal=FALSE,onefile=FALSE,paper="special")
postscript(file.name,width=6.4,height=7.0,horizontal=FALSE,onefile=FALSE,paper="special")

xrange <- c(-1.0,1.0,5)
xrange.axis <- c(-1.0,1.0,5)
yrange <- c(-1.0,1.0,5)
yrange.axis <- c(-1.0,1.0,5)

par(oma = c(0,0,0,0) )
source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange.R")
source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange_wy.R")
source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange_wox_wy.R")
source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange_woxyaxis.R")

source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap_wrange.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap_wrange_wy.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap_wrange_woxwy.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap_wrange_woxyaxis.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMapwrangewxwy_wpath.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMapwrangewxwoy_wpath.R")

source("~/papers/CG-FG_TACCM_REMD/FillConBar.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConBar.R")

nf <- layout(matrix(c(1,2,3,4,21,5,6,7,8,21,9,10,11,12,21,13,14,15,16,21,17,18,19,20,21),5,5,byrow=TRUE),c(4.5,3,3,3,1.5),c(1,1,1,1,1.4))

#par(mar =c(5,5,2,2) )
par(cex.axis=1.44)
par(cex.lab=1.44)
#label.x <- expression(paste(phi))
#label.y <- expression(paste(psi))

label.x <- " " # expression(paste(phi,"(radian/",pi,")"))
label.y <- expression(paste(psi,"(radian/",pi,")"))

name.pmf <- paste("${filename[1]}",sep='')
name.path <- paste("${filenamepath}",sep='')

xrange <- c(-1.0,1.0,4)
xrange.axis <- c(-1.0,1.0,4)
yrange <- c(-1.0,1.0,4)
yrange.axis <- c(-1.0,1.0,4)

cat("${filename[1]}",'\n')
felwFillConMapwrange(name.pmf,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",mar=c(0,5,1,0),xrange=xrange,yrange=yrange,plot.axis="yes",norm="T",xaflag="F",yaflag="T")

xrange <- c(-0.5,1.0,3)
xrange.axis <- c(-0.5,1.0,3)
yrange <- c(-1.0,1.0,4)
yrange.axis <- c(-1.0,1.0,5)

cat("${filename[2]}",'\n')
felwFillConMapwrange("${filename[2]}",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",mar=c(0,0,1,0),xrange=xrange,yrange=yrange,plot.axis="no",norm="T",xaflag="F")

cat("${filename[3]}",'\n')
felwFillConMapwrange("${filename[3]}",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",mar=c(0,0,1,0),xrange=xrange,yrange=yrange,plot.axis="no",norm="T",xaflag="F")

cat("${filename[4]}",'\n')
felwFillConMapwrange("${filename[4]}",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",mar=c(0,0,1,0),xrange=xrange,yrange=yrange,plot.axis="no",norm="T",xaflag="F")

cat("${filename[5]}",'\n')
felwFillConMapwrange("${filename[5]}",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",mar=c(0,5,0,0),xrange=xrange,yrange=yrange,plot.axis="yes",norm="T",xaflag="F",yaflag="T")

cat("${filename[6]}",'\n')
felwFillConMapwrange("${filename[6]}",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",mar=c(0,0,0,0),xrange=xrange,yrange=yrange,plot.axis="no",norm="T",xaflag="F")

cat("${filename[7]}",'\n')
felwFillConMapwrange("${filename[7]}",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",mar=c(0,0,0,0),xrange=xrange,yrange=yrange,plot.axis="no",norm="T",xaflag="F")

cat("${filename[8]}",'\n')
felwFillConMapwrange("${filename[8]}",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",mar=c(0,0,0,0),xrange=xrange,yrange=yrange,plot.axis="no",norm="T",xaflag="F")

cat("${filename[9]}",'\n')
felwFillConMapwrange("${filename[9]}",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",mar=c(0,5,0,0),xrange=xrange,yrange=yrange,plot.axis="yes",norm="T",xaflag="F",yaflag="T")

cat("${filename[10]}",'\n')
felwFillConMapwrange("${filename[10]}",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",mar=c(0,0,0,0),xrange=xrange,yrange=yrange,plot.axis="no",norm="T",xaflag="F")

cat("${filename[11]}",'\n')
felwFillConMapwrange("${filename[11]}",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",mar=c(0,0,0,0),xrange=xrange,yrange=yrange,plot.axis="no",norm="T",xaflag="F")

cat("${filename[12]}",'\n')
felwFillConMapwrange("${filename[12]}",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",mar=c(0,0,0,0),xrange=xrange,yrange=yrange,plot.axis="no",norm="T",xaflag="F")

cat("${filename[13]}",'\n')
felwFillConMapwrange("${filename[13]}",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",mar=c(0,5,0,0),xrange=xrange,yrange=yrange,plot.axis="yes",norm="T",xaflag="F",yaflag="T")

cat("${filename[14]}",'\n')
felwFillConMapwrange("${filename[14]}",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",mar=c(0,0,0,0),xrange=xrange,yrange=yrange,plot.axis="no",norm="T",xaflag="F")

cat("${filename[15]}",'\n')
felwFillConMapwrange("${filename[15]}",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",mar=c(0,0,0,0),xrange=xrange,yrange=yrange,plot.axis="no",norm="T",xaflag="F")

cat("${filename[16]}",'\n')
felwFillConMapwrange("${filename[16]}",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",mar=c(0,0,0,0),xrange=xrange,yrange=yrange,plot.axis="no",norm="T",xaflag="F")

xrange <- c(-1.0,0.5,3)
xrange.axis <- c(-1.0,0.5,3)
yrange <- c(-1.0,0.5,3)
yrange.axis <- c(-1.0,0.5,3)

cat("${filename[17]}",'\n')
felwFillConMapwrange("${filename[17]}",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",mar=c(5,5,0,0),xrange=xrange,yrange=yrange,plot.axis="yes",norm="T",xaflag="T",yaflag="T")

xrange <- c(-1.0,0.5,3)
xrange.axis <- c(-1.0,0.5,3)

cat("${filename[18]}",'\n')
felwFillConMapwrange("${filename[18]}",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",mar=c(5,0,0,0),xrange=xrange,yrange=yrange,plot.axis="no",norm="T")

xrange <- c(-1.0,0.5,3)
xrange.axis <- c(-1.0,0.5,3)

cat("${filename[19]}",'\n')
felwFillConMapwrange("${filename[19]}",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",mar=c(5,0,0,0),xrange=xrange,yrange=yrange,plot.axis="no",norm="T")

xrange <- c(-1.0,1.0,4)
xrange.axis <- c(-1.0,1.0,4)

cat("${filename[20]}",'\n')
felwFillConMapwrange("${filename[20]}",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",mar=c(5,0,0,0),xrange=xrange,yrange=yrange,plot.axis="no",norm="T")

par(mar=c(5,2,1,1))

felwFillConBar("${filename[1]}",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

EOF

Rscript ${parafile}

