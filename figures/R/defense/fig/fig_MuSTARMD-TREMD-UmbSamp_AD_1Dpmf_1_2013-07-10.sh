#!~/bin/sh

parafile1=~/calspa/TACCM_CGAAREMD/AD/para/2013-05-04/para_e_CG-FG_NR=8_2_TZ=750_fq=10ps_99SB_KZmax=5000_weljd0.001_mZ=100.sh

parafile2=~/calspa/refcalc/REMD/AD/para/para_s_REMD_vac_ff99SB_2012-07-24.sh
parafileUmbSam=~/calspa/refcalc/UmbSam/AD/para/para_s_vac_ff99SB_K=10_2012-11-12.sh
parafile_TAMD=~/calspa/TACCM/AD/para/para_e_ff99SB_NH_2012-07-31_TZ=750.sh
height=20 
width1=1 
numuene=4 

dir=~/calspa/TACCM_CGAAREMD/AD
dirOUT=~/defense

source ${parafile1}

height=${height} 

TZ=${TZs[1]}

AACGflag=CG

paraRfil=${dirOUT}/graph/fig_MuSTARMD-TREMD-UmbSamp_AD_1Dpmf_2013-07-10.R

cat <<EOF > ${paraRfil}
fact.y=1

TAA<-"${TAA}"
TCG<-"${TCG}"

TZs <- c( "${TZ}" )
nTZs <- length(TZs)

tau<-c(  "1.0" )
ntau<-length(tau)

mZ<-c(  "${mZ[1]}" )
nmZ<-length(mZ)

nKZAA<-${nKZAA}
nKZCG<-${nKZCG}
numRE<-${numRE}

pname1<-"${pname}"

numEX1<-"${numEX}"

TLbase1<-"${TLbase}"

AACG <- "${AACGflag}"
width <- "${width1}"

KZAAo <- "${KZAAo[2]}"
KZCGo <- "${KZCGo[2]}"

numuene <- "${numuene}"

height <- ${height}

phsiflag <- "phi"
phsi <- "phi"
name.title <- paste(AACG,sep="")

EOF

source ${parafile2}

cat <<EOF >> ${paraRfil}

ff<-"${ff}"

pname2<-"${pname}"

numEX2<-"${numEX}"

TLbase2<-"${TLbase}"
EOF

source ${parafileUmbSam}
TLns=`expr ${TLbase} / 1000`

cat <<EOF >> ${paraRfil}
ff2<-"${ff}"

pname3<-"${pname}"

TLns<-"${TLns}"

EOF

source ${parafile_TAMD}
cat <<EOF >> ${paraRfil}
ffTAMD <-"${ff}"
tauTAMD <-"1.0"
TBTAMD <- ${TB[1]}
KZTAMD <- "1000"
mZTAMD <- "50000.00"
TeTAMD=${T}

EOF

cat <<EOF >> ${paraRfil}

name.out <- "${dirOUT}/eps/fig_MuSTARMD-TREMD-UmbSamp_AD_1Dpmf_1_2013-07-10"

source("~/Rspa/plmGeneral_wsrange.R")

dir01 <- "/home/yamamori/calspa/TACCM_CGAAREMD/AD"
dirbase1 <- paste(dir01,"/e_CG-FG_NH_2012-07-27",sep="")

dir02 <- "~/calspa/refcalc/REMD/AD"
dirbase2 <- paste(dir02,"/s_REVAC_2012-07-23_",ff,sep="")

dir03 <- "~/calspa/refcalc/UmbSam/AD"
dirbase3 <- paste(dir03,"/s_UmbSam_vac_2012-11-12_",ff2,sep="")

dir04 <- "~/calspa/TACCM/AD"
dirbase4 <- paste(dir04,"/e_TACCM_NH_2012-07-31_",ffTAMD,sep="")

nTZs<-1

title=name.title

label.y="FEL (kT)"

xrange <- c(-180,180,6)
xrange.axis <- c(-180,180,6)

yrange <- c(0,12,6)     
yrange.axis <- c(0,12,6)

fact.x <- 1/pi*180.0
fact.y <- 1
ave.x <- 1

file.name <- paste(name.out,".eps",sep='')
#cat(file.name,'\n')
#tiff(file.name,width=800,height=700)
#tiff(file.name,width=600,height=525)
postscript(file.name,width=3.5,height=2.0,horizontal=FALSE,onefile=FALSE,paper="special")

par(mfrow=c(1,1))
par(mai=c(0.0,0.0,0.0,0.0))
par(oma = c(3.5,5.0,1.0,1.0) )
#par(oma=c(7.5,6.0,2.0,2.0))

s <- 1
name <- NULL

id.xs <- c(1,1,1,1)
id.ys <- c(2,2,2,2)
ids.ys <- c(3,3,3,3)
iro <- c(1,3,2,4)
senshu <- c(1,1,1,1)
#  tenshu <- c(1,1)
tenshu <- c(19,3,3,3)
hutosa <- c(0.8,0.8,0.8,0.8)
is.leg <- 0

phsiflag="phi"
xy <- "x"

#name[1] <- paste(dirbase1,"/tau=",tau[1],"/mZ=",mZ[1],"/TZ=",TZs[1],"/",pname1,"/nEX=",numEX1,"/fq=",TLbase1,"/pmf/1D/pmf_pymbar_1_TAA=",TAA,"_TCG=",TCG,"_TZ=",TZs[1],"_KZAAo=",KZAAo,"_KZCGo=",KZCGo,"_",AACG,phsiflag,"@",6,sep="")

name[1] <- paste(dirbase3,"/",pname3,"/pmf/1D/pmf_1D_UmbMD_vac_",TLns,"ns_",xy,"=",2,".txt",sep="")
name[2] <- paste(dirbase2,"/",pname2,"/nEX=",numEX2,"/freq=",TLbase2,"ps","/pmf/1D/pmf_p1",phsiflag,"@",2,sep="")
name[3] <- paste(dirbase1,"/tau=",tau[1],"/mZ=",mZ[1],"/TZ=",TZs[1],"/",pname1,"/nEX=",numEX1,"/fq=",TLbase1,"/pmf/1D_2013_05_12/pmf_pymbar_1_TAA=",TAA,"_TCG=",TCG,"_TZ=",TZs[1],"_KZAAo=",KZAAo,"_KZCGo=",KZCGo,"_",AACG,"_",xy,"=",2,".txt",sep="")
name[4] <- paste(dirbase4,"/tau=",tauTAMD,"/TB=",TBTAMD,"/KZ=",KZTAMD,"/mZ=",mZTAMD,"/pmf_1D/pmf_1D_T=",TeTAMD,"_x=",2,".txt",sep="")

label.x<-expression(paste(psi,"(radian/",pi,")"))

plmGeneralwsrange(data.names=name,
                  sd.names=name,
                  id.ys=id.ys,
                  ids.ys=ids.ys,
                  is.sen=rep(0,20),
                  label.size=0.5,axis.size=1.0,
                  iro=iro,axis.ft="F",is.header="T",
                  sdiro=iro,
                  xrange=xrange,yrange=yrange,
                  sdyrange=yrange,
                  warrow="T")
box(lwd=1.0)

axis(2,yaxp=yrange.axis,lwd=1.0,cex.axis=1.0)
mtext(outer=T,label.y,side=2,line=3.0,cex=1.0)

axis(1,xaxp=xrange.axis,lwd=1.0,cex.axis=1.0)
mtext(label.x,side=1,line=2.0,cex=1.0)

dev.off()
EOF

Rscript ${paraRfil}; echo ${paraRfil}
