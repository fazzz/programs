#!~/bin/sh

mode=2
#tiffwidth=420
#tiffheight=200
#tiffwidth=480
#tiffheight=200
tiffwidth=720 
tiffheight=300
height=20
width=2

width_pmf=0.3

proname=AD

source ~/calspa/TACCM_CGAAREMD/AD/para/2012-11-02/para_e_CG-FG_NR=4_TZ=700_fq=10ps_99SB_KZmax=5000_weljd0.001_mZ=100.sh

dir=~/calspa/MFEP/AD/

filenamepmf=~/calspa/TACCM_CGAAREMD/AD/e_CG-FG_NH_2012-07-27/tau=${tau[1]}/mZ=${mZ[1]}/TZ=${TZs[1]}/${pname}/nEX=${numEX}/fq=${TLbase}/pmf/pmf_TAA=${TAA}_TCG_${TCG}_TZ_${TZs[1]}_${width_pmf}_${KZAAo[$mode]}_${KZCGo[$mode]}_${MODES[$mode]}_4_2012-08-21

source ~/calspa/refcalc/REMD/AD/para/para_s_REMD_vac_ff99SB_2012-07-24.sh

filenamepmf2=~/calspa/refcalc/REMD/AD/s_REVAC_2012-07-23_${ff}/${pname}/nEX=${numEX}/freq=${TLbase}ps/pmf/pmf_ADv

parafile=~/gakkai/MSSJ_2012/graph/fig_poster_MuSTAMDcompTREMD2DFEL_2012-11-19.R
cat <<EOF > ${parafile}
source("~/calspa/GOLMAA/SH3/R/filled.contour.mod.R")
source("~/calspa/GOLMAA/SH3/R/fel3.R")
source("~/Rspa/plRamachad3.R")

dir <- "~/gakkai/MSSJ_2012"

title <- NULL

label.x <- expression(paste(phi))
label.y <- expression(paste(psi))

name.path <- paste("${filenamepath}",sep='')
name.out <- paste(dir,"/tiff/","fig_poster_MuSTAMDcompTREMD2DFEL_2012-11-19",sep='')
   
level <- seq(0.0,${height},${width})

file.name <- paste(name.out,'.tiff',sep='')
tiff(file.name,width=${tiffwidth},height=${tiffheight})

par(oma = c(0,0,0,0) )
source("~/papers/CG-FG_TACCM_REMD/FillConMap.R")
source("~/papers/CG-FG_TACCM_REMD/FillConBar.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConBar.R")

nf <- layout(matrix(c(1,2,3),1,3,byrow=TRUE),c(3,3,1),c(1))

#par(mar =c(5,5,2,2) )
#par(cex.axis=2.0)
#par(cex.lab=2.0)

par(mar =c(4,4,2,2) )
par(cex.axis=2.0)
par(cex.lab=2.0)

felwFillConMap("${filenamepmf2}",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

name.pmf <- paste("${filenamepmf}",sep='')
felwFillConMap(name.pmf,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

felwFillConBar(name.pmf,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

dev.off()
EOF

Rscript ${parafile}
