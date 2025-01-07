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

source ~/calspa/MFEP/AD/para/para_e_CG-FG_NR=4_wrefd0.1_TZ=750_fq=10ps_99SB.sh

dir=~/calspa/MFEP/AD/

filenamepmf=~/calspa/TACCM_CGAAREMD/AD/e_CG-FG_NH_2012-07-27/tau=${tau[1]}/mZ=${mZ[1]}/TZ=${TZs[1]}/${pname}/nEX=${numEX}/fq=${TLbase}/pmf/pmf_TAA=${TAA}_TCG_${TCG}_TZ_${TZs[1]}_0.3_${KZAAo[$mode]}_${KZCGo[$mode]}_${MODES[$mode]}_4_2012-08-21

filenameMGaussian=${dir}/pmf2MGaussian_2012-10-10_bo_e_CG-FG_NH_2012-07-27/tau=${tau[1]}/mZ=${mZ[1]}/TZ=${TZs[1]}/${pname}/nEX=${numEX}/fq=${TLbase}/pmf2MGaussian_TAA=${TAA}_TCG=${TCG}_${KZAAo[$mode]}_${KZCGo[$mode]}_${MODES[$mode]}_K=${K}.txt

filenamepath=${dir}/pmf2MGaussian_2012-10-10_bo_e_CG-FG_NH_2012-07-27/tau=${tau[1]}/mZ=${mZ[1]}/TZ=${TZs[1]}/${pname}/nEX=${numEX}/fq=${TLbase}/path_TAA=${TAA}_TCG=${TCG}_${KZAAo[$mode]}_${KZCGo[$mode]}_${MODES[$mode]}_K=${K}_@${xi},${yi}-@${yi},${yf}.txt

filenamepath2=${dir}/pmf2MGaussian_2012-10-10_bo_e_CG-FG_NH_2012-07-27/tau=${tau[1]}/mZ=${mZ[1]}/TZ=${TZs[1]}/${pname}/nEX=${numEX}/fq=${TLbase}/path2_TAA=${TAA}_TCG=${TCG}_${KZAAo[$mode]}_${KZCGo[$mode]}_${MODES[$mode]}_K=${K}_@${xi},${yi}-@${yi},${yf}.txt

parafile=~/papers/Path_Search/graph/fig_1Dpmf_along_MFEP_fMuSTARMD_2012-10-24.R
cat <<EOF > ${parafile}
source("~/calspa/GOLMAA/SH3/R/filled.contour.mod.R")
source("~/calspa/GOLMAA/SH3/R/fel3.R")
source("~/Rspa/plRamachad3.R")

dir <- "~/papers/Path_Search"

title <- NULL

label.x <- expression(paste(phi))
label.y <- expression(paste(psi))

name.path <- paste("${filenamepath}",sep='')
name.out <- paste(dir,"/tiff/","fig_1Dpmf_along_MFEP_fMuSTARMD_2012-10-24",sep='')
   
level <- seq(0.0,${height},${width})

file.name <- paste(name.out,'.tiff',sep='')
tiff(file.name,width=${tiffwidth},height=${tiffheight})

par(oma = c(0,0,0,0) )
source("~/papers/CG-FG_TACCM_REMD/FillConMap.R")
source("~/papers/CG-FG_TACCM_REMD/FillConBar.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap_wpath.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConBar.R")

nf <- layout(matrix(c(1,1,2,3),2,2,byrow=TRUE),c(3,1),c(1,1))

par(mar =c(5,5,2,2) )
par(cex.axis=2.0)
par(cex.lab=2.0)

label.x=" "
label.y="pmf"

xrange <- c(0,1.0,10)
xrange.axis <- xrange
yrange <- c(0,${height},10)
yrange.axis <- yrange

fact.x <- 1
fact.y <- 1
ave.x <- 1

id.xs <- c(1,1)
id.ys <- c(2,2)
ids.ys <- c(3,3)
iro <- c(1,2)
senshu <- c(1,1)
tenshu <- c(3,3)
hutosa <- c(1,1)
is.leg <- 0

source("~/Rspa/plmGeneral_wsrange.R")

########################################################
# plmGeneralwsrange(data.names="${filenamepath2}",     #
#                   sd.names="${filenamepath2}",       #
#                   id.ys=id.ys,		       #
#                   ids.ys=ids.ys,		       #
#                   is.sen=rep(0,20),		       #
#                   label.size=0.5,axis.size=2.0,      #
#                   iro=iro,axis.ft="F",is.header="T", #
#                   sdiro=iro,			       #
#                   xrange=xrange,yrange=yrange,       #
#                   sdyrange=yrange,		       #
#                   warrow="T")			       #
########################################################
#box(lwd=2.0)
  
#axis(2,yaxp=yrange.axis,lwd=2.0,cex.axis=1.5)
#mtext(outer=T,label.y,side=2,line=4.0,cex=1.5)
  
#axis(1,xaxp=xrange.axis,lwd=2.0,cex.axis=1.5)
#mtext(outer=T,label.x,side=1,line=4.0,cex=1.5)

name.pmf <- paste("${filenamepmf}",sep='')
cat(name.pmf," ",name.path,'\n')
felwFillConMapwpath(name.pmf,name.path,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

name.pmf <- paste("${filenamepmf}",sep='')
cat(name.pmf," ",name.path,'\n')
felwFillConMapwpath(name.pmf,name.path,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

felwFillConBar(name.pmf,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

dev.off()
EOF

Rscript ${parafile}
