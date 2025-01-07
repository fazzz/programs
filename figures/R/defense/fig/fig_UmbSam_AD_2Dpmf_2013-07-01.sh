#!~/bin/sh

proname=AD
dirUMB=~/calspa/refcalc/UmbSam/AD
dirOUT=~/defense/

parafile=~/calspa/refcalc/UmbSam/AD/para/para_s_vac_ff99SB_K=10_2012-11-12.sh
source ${parafile}

nb=21

TLns=`expr ${TLbase} / 1000`
direqubase=${dirUMB}/s_UmbSam_vac_2012-11-12_${ff}/
dirpmf=${direqubase}${pname}/pmf
filenamepmf=${dirpmf}/pmf_UmbMD_vac_${TLns}ns.txt

paraRfile=${dirOUT}/graph/para_UmbSam_AD_2Dpmf_2013-07-01.R

cat <<EOF > ${paraRfile}
fact.x <- 180/pi
fact.y <- 180/pi
#fact.p <- 1.0/pi

dir    <- "${dirUMB}"
dirout <- "${dirOUT}"

title <- NULL

label.x <- expression(paste(""))
label.y <- expression(paste(""))

name.out <- paste(dirout,"/eps/","fig_UmbSam_AD_2Dpmf_2013-07-01",sep='')
   
level <- seq(0.0,20.0,2.0)

xrange <- c(-180,180,6)
xrange.axis <- c(-180,180,6)
yrange <- c(-180,180,6)
yrange.axis <- c(-180,180,6)

file.name <- paste(name.out,'.eps',sep='')
postscript(file.name,width=4.1496,height=2.6472,horizontal=FALSE,onefile=FALSE,paper="special")

par(oma = c(0,0,0,0) )
source("~/defense/fig/FillConMap.R")
source("~/defense/fig/FillConBar.R")
source("~/defense/fig/fel_FillConMap.R")
source("~/defense/fig/fel_FillConBar.R")

nf <- layout(matrix(c(1,2),1,2,byrow=TRUE),c(6.47,3.07),c(1))

par(cex.axis=1.2)
par(cex.lab=1.2)

name <- paste("${filenamepmf}_2",sep='')

cat(name)

label.x <- ""
label.y <- ""

felwFillConMap(name,label.x=label.x,label.y=label.y,level=level,
              title=title,kt="F/kBT",
              xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,
              mai=c(0.5,0.5,0.1,0.0))

felwFillConBar(name,label.x=label.x,label.y=label.y,level=level,
                   mai=c(1.348,0.75,0.1,0.75))

EOF

Rscript ${paraRfile}; echo ${paraRfile}
