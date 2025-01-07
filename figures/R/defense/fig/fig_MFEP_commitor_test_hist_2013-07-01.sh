#!~/bin/sh

K=19
mode=2

xi=-1.4
yi=1.1
xf=1.1
yf=-0.8

minx=-2.5
maxx=2.5
miny=-2.5
maxy=2.5

source ~/calspa/MFEP/AD/para/para_e_CG-FG_NR=4_wrefd0.1_TZ=750_fq=10ps_99SB.sh

crdTSL=ADv_TS_-0.019--1.010.crd_3_1

pointsize=0.5

dir=~/calspa/MFEP/AD
dirOUT=~/defense

filename=( dummy )

filename=${dir}/committor_test_2013-03-24_bo_e_CG-FG_NH_2012-07-27/tau_${tau[1]}/mZ_${mZ[1]}/TZ_${TZs[1]}/${pname}/nEX_${numEX}/fq_${TLbase}/${crdTSL}/anl/comm_hist.txt 

paraRfil=${dirOUT}/graph/graph_fig_6_b_committor_test_hist_2013-04-10.R

cat <<EOF > ${paraRfil}
setwd("~/Rspa")
source("~/Rspa/ini.R")
source("~/Rspa/plmGeneral_wsrange_wdrange_wdiffx.R")

name.title <- paste("${dirOUT}/eps/fig_MFEP_commitor_test_hist_2013-07-01",sep="")

leg.pos="topright"
is.leg <- 0

xrange <- c(0,102,102)
xrange.axis <- c(0,100,10)
yrange <- c(0.0,50,5)
yrange.axis <- c(0.0,50,5)

file.name <- paste(name.title,".eps",sep="")
postscript(file.name,width=4.2,height=3.4,horizontal=FALSE,onefile=FALSE,paper="special")

par(mfrow=c(1,1))
par(mar=c(0.0,0.0,0.0,0.0))
par(oma=c(7.5,6.0,2.0,2.0))

label.x <- expression(paste("index of image of string"))
label.y <- "number of histgram"

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
                               width=10.0,
                               seihouflag=T,
                               warrow="F",ave.x=ave.x)

box(lwd=2.0)
          
axis(1,xaxp=xrange.axis,lwd=2.0,cex=1.0)
mtext(outer=T,label.x,side=1,line=3.0,cex=1.0)
        
axis(2,yaxp=yrange.axis,lwd=2.0,cex=0.5)
mtext(label.y,side=2,line=2.0,cex=1.0)

EOF

Rscript ${paraRfil}

