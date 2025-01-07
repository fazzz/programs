#!~/bin/sh

mode=2

minx=-4.0
maxx=4.0
miny=-4.0
maxy=4.0

source ~/calspa/MFEP/AD/para/para_e_CG-FG_NR=8_2_TZ=750_fq=10ps_99SB_KZmax=5000_weljd0.001_mZ=100_2013-06-05.sh

pointsize=0.5

dir=~/calspa/MFEP/AD
dirOUT=~/papers/Path_Search

filename=( dummy )

filename=${dir}/pmf2MGaussian_2012-10-10_bo_e_CG-FG_NH_2013-05-04/tau=${tau[1]}/mZ=${mZ[1]}/TZ=${TZs[1]}/${pname}/nEX=${numEX}/fq=${TLbase}/KLdv_TAA=${TAA}_TCG=${TCG}_${KZAAo[$mode]}_${KZCGo[$mode]}_${MODES[$mode]}_wx${minx}-${maxx}_wy${miny}-${maxy}.txt

paraRfil=${dirOUT}/graph/graph_fig_4_KLdiv_K_2013-06-14.R

cat <<EOF > ${paraRfil}
setwd("~/Rspa")
source("~/Rspa/ini.R")
source("~/Rspa/plmGeneral_wsrange_wdrange_wdiffx.R")

name.title <- paste("${dirOUT}/eps/fig_4_KLdiv_K_2013-06-14",sep="")

leg.pos="topright"
is.leg <- 0

xrange <- c(0,20,10)
xrange.axis <- c(0,20,10)
yrange <- c(0.0,2.1,5)
yrange.axis <- c(0.0,2.0,4)

file.name <- paste(name.title,".eps",sep="")
postscript(file.name,width=3.5,height=2.8,horizontal=FALSE,onefile=FALSE,paper="special")
#postscript(file.name,width=5.2,height=4.2,horizontal=FALSE,onefile=FALSE,paper="special")

par(mfrow=c(1,1))
par(mar=c(0.0,0.0,0.0,0.0))
par(oma=c(7.5,6.0,2.0,2.0))

label.x <- expression(paste("K "))
label.y <- "KL divergence"

sen<-c(0)
senshu <- 1
iro<-c(2)
id.ys<-2

filex.names<-paste("${filename}",sep="")
filey.names<-paste("${filename}",sep="")

hutosa <- 1.0
ave.x <- 1
          
plmGeneralwsrangewdrangewdiffx(datax.names=filex.names,
                               datay.names=filey.names,
                               sd.names=file.names,
                               id.ys=id.ys,
                               ids.ys=ids.ys,
                               label.size=0.5,axis.size=2.0,
                               iro=iro,axis.ft="F",is.header="F",
                               sdiro=iro,
                               xrange=xrange,yrange=yrange,
                               sdyrange=yrange,
                               is.sen=sen,width=10.0,
                               warrow="F",ave.x=ave.x,
                               point.size=${pointsize})

box(lwd=2.0)
          
axis(1,xaxp=xrange.axis,lwd=2.0,cex=1.0)
mtext(outer=T,label.x,side=1,line=3.0,cex=1.0)
        
axis(2,yaxp=yrange.axis,lwd=2.0,cex=0.5)
mtext(label.y,side=2,line=2.0,cex=1.0)

EOF

Rscript ${paraRfil}
