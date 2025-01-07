setwd("~/Rspa")
source("~/Rspa/ini.R")
source("~/Rspa/plmGeneral_wsrange_wdrange_wdiffx.R")

name.title <- paste("/home/yamamori/papers/Path_Search/eps/fig_6_MFEP1D-Met-Enkephalin_2014-04-10",sep="")

leg.pos="topright"
is.leg <- 0

xrange <- c(0,1.0,10)
xrange.axis <- c(0,1.0,10)
yrange <- c(0.0,4,4)
yrange.axis <- c(0.0,4,4)

file.name <- paste(name.title,".eps",sep="")
#postscript(file.name,width=3.5,height=2.8,horizontal=FALSE,onefile=FALSE,paper="special")
#postscript(file.name,width=5.2,height=4.2,horizontal=FALSE,onefile=FALSE,paper="special")
postscript(file.name,width=4.2,height=3.4,horizontal=FALSE,onefile=FALSE,paper="special")

par(mfrow=c(1,1))
par(mar=c(0.0,0.0,0.0,0.0))
par(oma=c(7.5,6.0,2.0,2.0))

label.x <- expression(paste(""))
label.y <- "MFEP(/KT)"

sen<-c(0)
senshu <- 1
iro<-c(2)
id.ys<-2

filex.names<-paste("/home/yamamori/calspa/MFEP/MetEnk/pmf2MGaussian_2013-04-10_bo_e_CG-FG_CMD_NH_2012-06-04/path2_TAA=300_TCG1=370_TCG2=370_1000_0_0_K=12_@0.06,0.15-@0.17,0.06.txt",sep="")
filey.names<-paste("/home/yamamori/calspa/MFEP/MetEnk/pmf2MGaussian_2013-04-10_bo_e_CG-FG_CMD_NH_2012-06-04/path2_TAA=300_TCG1=370_TCG2=370_1000_0_0_K=12_@0.06,0.15-@0.17,0.06.txt",sep="")

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
                               point.size=0.5)

box(lwd=2.0)
          
axis(1,xaxp=xrange.axis,lwd=2.0,cex=1.0)
mtext(outer=T,label.x,side=1,line=3.0,cex=1.0)
        
axis(2,yaxp=yrange.axis,lwd=2.0,cex=0.5)
mtext(label.y,side=2,line=2.0,cex=1.0)

