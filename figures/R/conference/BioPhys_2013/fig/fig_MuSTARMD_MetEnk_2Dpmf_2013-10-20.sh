#!~/bin/sh

proname=MetEnk
dirOUT=~/gakkai/BioPhys_2013

paraRfile=${dirOUT}/graph/para_MuSTARMD_MetEnk_2Dpmf_2013-10-20.R

cat <<EOF > ${paraRfile}

name <- "~/calspa/TACCM_CGAAREMD/MetEnk/e_CG-FG_CMD_NH_2013-08-31/1FG2CG/tau=1.0/mZ=1000.00/ep=0.01/cutoff=4.7/TZ=1000/TCG=400/2CG_2013-09-05_KAA=1250_KCG=1000-500_2/freq=1/fb=1.0/fa=1.0/ft=0.2/fc=1.0/fn=1.0/pmf/pmf_TAA=300_TCG1_400_TCG2_400_TZ_1000_0.015_1250_0_0_AA_bo10000ps_2013-08-24_2013-08-31"

fact.x <- 1
fact.y <- 1

dirout <- "${dirOUT}"

title <- NULL

label.x <- expression(paste(""))
label.y <- expression(paste(""))

name.out <- paste(dirout,"/eps/","fig_MuSTARMD_MetEnk_2Dpmf_2013-10-20",sep='')
   
level <- seq(0,20,2)

xrange <- c(0,0.3,3)       
xrange.axis <- c(0,0.3,3)  
yrange <- c(0,0.3,3)       
yrange.axis <- c(0.1,0.3,2)

file.name <- paste(name.out,'.eps',sep='')
postscript(file.name,width=4.0,height=2.5512,horizontal=FALSE,onefile=FALSE,paper="special")

par(oma = c(0,0,0,0) )
source("~/defense/fig/FillConMap.R")
source("~/defense/fig/FillConBar.R")
source("~/defense/fig/fel_FillConMap.R")
source("~/defense/fig/fel_FillConBar.R")

nf <- layout(matrix(c(1,2),1,2,byrow=TRUE),c(6.47,3.07),c(1))

par(cex.axis=1.2)
par(cex.lab=1.2)

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
