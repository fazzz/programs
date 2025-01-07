#!~/bin/sh

proname=AD
dirUMB=~/calspa/refcalc/UmbSam/AD
dirOUT=~/papers/CG-FG_TACCM_REMD

parafile=~/calspa/refcalc/UmbSam/AD/para/para_s_vac_ff99SB_K=10_2012-11-12.sh
source ${parafile}

tiffwidth=400
tiffheight=300

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

name.out <- paste(dirout,"/tiff/","fig6_UmbSamp_AD_2Dpmf_2013-01-14",sep='')
   
level <- seq(0.0,10.0,1)

file.name <- paste(name.out,'.tiff',sep='')
tiff(file.name,width=${tiffwidth},height=${tiffheight})

par(oma = c(0,0,0,0) )
source("~/papers/CG-FG_TACCM_REMD/FillConMap.R")
source("~/papers/CG-FG_TACCM_REMD/FillConBar.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConBar.R")

nf <- layout(matrix(c(1,2),1,2,byrow=TRUE),c(3,1),c(1,1))

par(mar =c(5,5,2,2) )
par(cex.axis=2.0)
par(cex.lab=2.0)

name.pmf <- paste("${filenamepmf}",sep='')

felwFillConMap(name.pmf,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")
felwFillConBar(name.pmf,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

dev.off()
EOF

Rscript ${paraRfile}; echo ${paraRfile}


