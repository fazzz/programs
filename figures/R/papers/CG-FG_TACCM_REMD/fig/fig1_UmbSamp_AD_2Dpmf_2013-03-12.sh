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
dir    <- "${dirUMB}"
dirout <- "${dirOUT}"

title <- NULL

label.x <- expression(paste(phi))
label.y <- expression(paste(psi))

name.out <- paste(dirout,"/eps/","fig1_UmbSamp_AD_2Dpmf_2013-02-28",sep='')
   
level <- seq(0.0,10.0,1.0)
#level <- seq(0.0,20.0,2)

xrange <- c(-3,3,6)
yrange <- c(-3,3,6)

file.name <- paste(name.out,'.eps',sep='')
#tiff(file.name,width=${tiffwidth},height=${tiffheight})
postscript(file.name,width=3.2,height=2.6,horizontal=FALSE,onefile=FALSE,paper="special")

par(oma = c(0,0,0,0) )
source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange.R")
source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange_wy.R")
source("~/papers/CG-FG_TACCM_REMD/FillConBar.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap_wrange.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap_wrange_wy.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConBar.R")

nf <- layout(matrix(c(1,2),1,2,byrow=TRUE),c(3.0,1.0),c(1))

par(cex.axis=0.8)
par(cex.lab=1.0)

name.pmf <- paste("${filenamepmf}",sep='')

cat(name.pmf,"\n")

felwFillConMapwrangewy(name.pmf,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,mar=c(4,4,1,0))

par(mar=c(4,2,1,1))

felwFillConBar(name.pmf,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

#dev.off()
EOF

Rscript ${paraRfile}; echo ${paraRfile}


