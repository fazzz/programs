#!~/bin/sh

proname=MetEnk
dirOUT=~/defense/

parafile=~/calspa/TACCM_CGAAREMD/MetEnk/para/2013-05-31/para_e_CG-FG_NH_2CG_KAA=1750_TZ=2000_freq=1ps_2013-05-31.sh

source ${parafile}

paraRfile=${dirOUT}/graph/para_MuSTARMD_MetEnk_2Dpmf_2013-07-01.R

cat <<EOF > ${paraRfile}
max <- "0.3"

TAA<-"${TAA}"
TCG1<-"${TCG1}"
TCG2<-"${TCG2}"

TZs <-"${TZs[1]}"
tau<-"${tau[1]}"

ep<-"${ep}"
nb<-"${nb}"
cutoff<-"${cutoff}"

mZ<-"${mZ[1]}"

nKZAA<-${nKZAA}
nKZCG1<-${nKZCG1}
nKZCG2<-${nKZCG2}
numRE<-${numRE}

pname<-"${pname}"

freq<-"${TLbase}"

AACG <- "AA"
width <- "0.01"

KZAAo <- "${KZAAo[1]}"
KZCG1o <- "${KZCG1o[1]}"
KZCG2o <- "${KZCG2o[1]}"

level <- seq(0,10,1)

state <- "1FG2CG"

TS <- "${numEX}"

TLbase <- "${TLbase}"

dir0 <- "/home/yamamori/calspa/TACCM_CGAAREMD/MetEnk/"
dirbase <- paste(dir0,"/e_CG-FG_CMD_NH_2012-06-04/",state,sep="")

nTZs<-1

name <- paste(dirbase,"/tau=",tau,"/mZ=",mZ,"/ep=",ep,"/cutoff=",cutoff,"/TZ=",TZs,"/",pname,"/freq=",freq,"/pmf/pmf_pymbar_AA_41_2",sep="")

fact.x <- 1
fact.y <- 1

dirout <- "${dirOUT}"

title <- NULL

label.x <- expression(paste(""))
label.y <- expression(paste(""))

name.out <- paste(dirout,"/eps/","fig_MuSTARMD_MetEnk_2Dpmf_2013-07-01",sep='')
   
level <- seq(0,10,1)

xrange <- c(0,0.3,3)       
xrange.axis <- c(0,0.3,3)  
yrange <- c(0,0.3,3)       
yrange.axis <- c(0.1,0.3,2)

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
