#!~/bin/sh

proname=AD
dirOUT=~/defense/

parafile=~/calspa/TACCM/AD/para/para_e_ff99SB_NH_2012-07-31_TZ=750.sh

source ${parafile}

paraRfile=${dirOUT}/graph/para_TAMD_AD_2Dpmf_2013-07-01.R

i=1
cat <<EOF > ${paraRfile}
T_TAMD<-"${T}"

TB_TAMD<-"${TB[$i]}"

tau_TAMD<-"${tau}"

KZ_TAMD<-"${KZ}"

mZ_TAMD<-"${mZ}"

width_TAMD <- "0.3"

ff_TAMD <- "${ff}"

dir0 <- "/home/yamamori/calspa/TACCM/AD"
dirbase <- paste(dir0,"/e_TACCM_NH_2012-07-31_",ff_TAMD,sep="")
name <- paste(dirbase,"/tau=",tau_TAMD,"/TB=",TB_TAMD,"/KZ=",KZ_TAMD,"/mZ=",mZ_TAMD,"/pmf/pmf_T=",T_TAMD,".Zhist",sep="")

fact.x <- 180/pi
fact.y <- 180/pi
#fact.p <- 1.0/pi

dirout <- "${dirOUT}"

title <- NULL

label.x <- expression(paste(""))
label.y <- expression(paste(""))

name.out <- paste(dirout,"/eps/","fig_TAMD_AD_2Dpmf_2013-07-01",sep='')
   
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
