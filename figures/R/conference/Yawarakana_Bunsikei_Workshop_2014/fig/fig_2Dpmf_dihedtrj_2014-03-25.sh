#!~/bin/sh

dirOUT=~/gakkai/Yawarakana_Bunsikei_Workshop_2014/

dirpmf=/home/yamamori/data_paper/MuSTAR_MD/fig_3
dirdtrj=/home/yamamori/data_paper/MuSTAR_MD/fig_2

name=( dummy MuSTARMD TAMD REMD CMD )

filepmf=( dummy ${dirpmf}/pmf_AD_MuSTAR_MD ${dirpmf}/pmf_AD_TAMD ${dirpmf}/pmf_AD_REMD ${dirpmf}/pmf_AD_CMD )

filedtrj=( dummy ${dirdtrj}/AD_MuSTARMD.dihed_trj ${dirdtrj}/AD_TAMD.dihed_trj ${dirdtrj}/AD_REMD.dihed_trj ${dirdtrj}/AD_CMD.dihed_trj )

num=$(expr ${#name[*]} - 1)

for  i in `seq 1 ${num}`; do
    paraRfile=${dirOUT}/graph/fig_2Dpmf_${name[$i]}_2014-03-25.R

    cat <<EOF > ${paraRfile}

fact.x <- 180/pi
fact.y <- 180/pi

dirout <- "${dirOUT}"

title <- NULL

label.x <- expression(paste(""))
label.y <- expression(paste(""))

name.out <- paste(dirout,"/eps/","fig_AD_2Dpmf_${name[$i]}_2014-03-25",sep='')
   
level <- seq(0.0,30.0,2.0)

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

nf <- layout(matrix(c(1,2),1,2,byrow=TRUE),c(3,1.2),c(1))

par(cex.axis=1.2)
par(cex.lab=1.2)

label.x <- ""
label.y <- ""

name<-paste("${filepmf[$i]}",sep='')
cat(name)
felwFillConMap(name,label.x=label.x,label.y=label.y,level=level,
              title=title,kt="F/kBT",
              xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,
              mai=c(0.5,0.5,0.1,0.0))

felwFillConBar(name,label.x=label.x,label.y=label.y,level=level,
               mai=c(1.348,0.75,0.1,0.75))

EOF

    Rscript ${paraRfile}; echo ${paraRfile}
done

for  i in `seq 1 ${num}`; do
    paraRfile=${dirOUT}/graph/fig_dihedtrj_${name[$i]}_2014-03-25.R

    cat -n ${filedtrj[$i]} > ${dirOUT}/AD_${nama[$i]}.dihed_trj_n

    cat <<EOF > ${paraRfile}
setwd("~/Rspa")
source("~/Rspa/ini.R")
source("~/Rspa/plmGeneral_wsrange_wdrange_wdiffx.R")

dirout <- "${dirOUT}"

name.out <- paste(dirout,"/eps/","fig_AD_dihedtrj_${name[$i]}_2014-03-25",sep='')

leg.pos="topright"
is.leg <- 0

fact.y=1/pi

file.name <- paste(name.out,".eps",sep="")
postscript(file.name,width=1.5843,height=2.3716,horizontal=FALSE,onefile=FALSE,paper="special")

par(mfrow=c(2,1))
par(mar=c(0.0,0.0,0.0,0.0))
par(oma=c(0.5,0.5,0.1,0.1))

par(cex=1.0)

filex.names<-NULL
filey.names<-NULL

xrange <- c(0,10000,5)
xrange.axis <- c(0,10000,5)

yrange <- c(-1.0,1.0,4)
yrange.axis <- c(-1.0,1.0,4)

label.x <- expression(paste("time (ps)"))
label.y <- expression(paste("Dihedral Angle"))

hutosa <- NULL

for (i in 1:2) {
    sen<-c(0)
    senshu <- 1
    iro<-c(2)
    id.ys<-i

    filex.names<-"${dirOUT}/AD_${nama[$i]}.dihed_trj_n"
    filey.names<-"${filedtrj[$i]}"

    hutosa <- 0.1
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
          
    if ( i == 2 ) {
        axis(1,xaxp=xrange.axis,lwd=2.0,cex=1.0)
        mtext(outer=T,label.x,side=1,line=3.0,cex=1.0)
    }
        
    axis(2,yaxp=yrange.axis,lwd=2.0)
    mtext(outer=T,label.y,side=2,line=3.0,cex=1.0)
}

EOF

    Rscript ${paraRfile}; echo ${paraRfile}
done
