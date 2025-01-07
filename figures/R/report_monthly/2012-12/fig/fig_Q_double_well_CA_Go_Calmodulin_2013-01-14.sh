#!~/bin/sh

opt=(dummy width height  )
nopt=${#opt[*]}
if [ $# -le `expr ${nopt} - 2`  ]; then
    echo "USAGE " $0  ${opt[*]:1:${nopt}}
    echo $*
    exit
fi

num=1
while [ $num -le `expr ${nopt} - 1` ]; do
    eval ${opt[$num]}=$1
    shift 1
    num=`expr $num + 1`
done

dir=~/calspa/CaM/Go-like_model/MB/
dirout=~/Report/2012-12/

parafile=~/calspa/CaM/Go-like_model/MB/para/para_e_vT_vpara_GOLMCA_MB_JMB2001_NH_2013-01-14.sh
source ${parafile}

kaimin=-20 
kaimax=20

paraRfil=${dirout}/graph/graph_Q_double_well_CA_Go_Calmodulin_2013-01-14.R

cat <<EOF > ${paraRfil}
proname <- "${proname}"

state <- c( "fapo", "fholo" )
nstate <- length(state)

T <- "${T[1]}"

tau <- "${tau[1]}"

cutoff <- "${cutoff[1]}"

ep <- "${ep[1]}"

de <- "${de[1]}"

d <- "${d[1]}"

kaimin=${kaimin}
kaimax=${kaimax}

name.title <- paste("${dirout}/tiff/fig_Q_double_well_CA_Go_Calmodulin_2013-01-14",sep="")

setwd("~/Rspa")

source("~/Rspa/ini.R")
source("~/Rspa/plmGeneral_wsrange_wdrange_wdiffx.R")

leg.pos="topright"
is.leg <- 0

xrange <- c(0,100000000,4)
xrange.axis <- c(10000000,100000000,3)

label.x <- expression(paste("time step "))

file.name <- paste(name.title,".tiff",sep="")
tiff(file.name,width=${width},height=${height})

par(mfrow=c(2,1))
par(mar=c(0.0,0.0,0.0,0.0))
par(oma=c(7.5,6.0,2.0,2.0))

filex.names<-NULL
filey.names<-NULL

yrange <- c(kaimin,kaimax,10)

kaimin2 <- kaimin*0.9
kaimax2 <- kaimax*0.9
yrange.axis <- c(kaimin2,kaimax2,8)

label.y <- expression(paste(kai))

id.ys<-NULL

dir <- "~/calspa/CaM/Go-like_model/MB/e_vT_vpara_GOLMCA_MB_JMB2001_NH_2012-07-23/"

hutosa <- NULL

for (i in 1:nstate) {
    sen<-c(1)
    senshu <- 1
    iro<-c(2)
    id.ys<-c(1)

    dirbase.name <- paste(dir,"/",state[i],"/tau=",tau,"/ep=",ep,"/coff=",cutoff,"/de=",de,"/d=",d,"/",sep="")
        
    filex.names<-paste(dirbase.name,proname,"_equ_T=",T,".out",sep='')
    filey.names<-paste(dirbase.name,"anl/",proname,"_equ_T=",T,".kai",sep='')

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
                                   warrow="F",ave.x=ave.x)
    box(lwd=2.0)
          
    if ( i == 2 ) {
        axis(1,xaxp=xrange.axis,lwd=2.0,cex=1.0)
        mtext(outer=T,label.x,side=1,line=5.0,cex=2.0)
    }
        
    axis(2,yaxp=yrange.axis,lwd=2.0,cex=1.0)
    mtext(label.y,side=2,line=3.0,cex=2.0)
}

dev.off()

EOF

Rscript ${paraRfil}
