#!~/bin/sh

parafileMuSTARMD=~/calspa/TACCM_CGAAREMD/AD/para/2013-08-28/para_e_CG-FG_NR=8_TZ=1000_fq=10ps_99SB_KZmax=1000_weljd0.001_mZ=100.sh
parafileTREMD=~/calspa/refcalc/REMD/AD/para/para_s_REMD_vac_ff99SB_2013-08-17_TZ=300-1200.sh
parafileUmbSam=~/calspa/refcalc/UmbSam/AD/para/para_s_vac_ff99SB_K=10_2012-11-12.sh
parafileREUS=~/calspa/refcalc/REUS/AD/para/2013-08-30/para_s_REUSMD_vac_ff99SB_Nbin=4x4_K=0.5_5ns_freq=10ps_2013-08-30.sh
parafile_TAMD=~/calspa/TACCM/AD/para/2013-08-20/para_e_ff99SB_2013-08-20_mZ=100_Tz=800.sh
height=20 
width1=1 
numuene=4 

dir=~/calspa/TACCM_CGAAREMD/AD
dirOUT=~/papers/CG-FG_TACCM_REMD

source ${parafileMuSTARMD}

height=${height} 

TZ=${TZs[1]}

AACGflag=CG

paraRfil=${dirOUT}/graph/para_fig4-a-b_comp_MuSTARMD-TREMD-REUS-TAMD-UmbSamp_AD_1Dpmf_2013-09-02.R

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

source ${parafileTREMD}

cat <<EOF >> ${paraRfil}

ff<-"${ff}"

pname2<-"${pname}"

numEX2<-"${numEX}"

TLbase2<-"${TLbase}"
EOF

source ${parafileREUS}

cat <<EOF >> ${paraRfil}
TLbaseREUS<-"${TLbase}"

numEXREUS<-"${numEX}"

numREUS<-"${numREUS}"

pnameREUS<-"${pname}"
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
mZTAMD <- "100.00"
TeTAMD=${T}

EOF

cat <<EOF >> ${paraRfil}

name.out <- "${dirOUT}/eps/fig4-a-b_comp_MuSTARMD-TREMD-REUS-TAMD-UmbSamp_AD_1Dpmf_2013-09-02"

source("~/Rspa/plmGeneral_wsrange.R")

dir01 <- "/home/yamamori/calspa/TACCM_CGAAREMD/AD"
dirbase1 <- paste(dir01,"/e_CG-FG_NH_2012-07-27",sep="")

dir02 <- "~/calspa/refcalc/REMD/AD"
dirbase2 <- paste(dir02,"/s_REVAC_2012-07-23_",ff,sep="")

dirREUS <- "~/calspa/refcalc/REUS/AD"
dirbaseREUS <- paste(dirREUS,"/s_REUSVAC_2013-08-30_",ff,sep="")

dir03 <- "~/calspa/refcalc/UmbSam/AD"
dirbase3 <- paste(dir03,"/s_UmbSam_vac_2012-11-12_",ff2,sep="")

dir04 <- "~/calspa/TACCM/AD"
dirbase4 <- paste(dir04,"/e_TACCM_NH_2012-07-31_",ffTAMD,sep="")

nTZs<-1

title=name.title

label.y="FEL (kT)"

xrange <- c(-1.0,1.0,4)
xrange.axis <- c(-1.0,1.0,4)

yrange <- c(-0.1,12,6)
yrange.axis <- c(0,12,6)

fact.x <- 1/pi
ave.x <- 1

file.name <- paste(name.out,".eps",sep='')
#cat(file.name,'\n')
#tiff(file.name,width=800,height=700)
#tiff(file.name,width=600,height=525)
postscript(file.name,width=2.7,height=3.0,horizontal=FALSE,onefile=FALSE,paper="special")

par(mfrow=c(2,1))
par(mai=c(0.0,0.0,0.0,0.0))
par(oma = c(3.5,3.0,1.0,0.6) )
#par(oma=c(7.5,6.0,2.0,2.0))

s <- 1
name <- NULL

id.xs <- c(1,1,1,1,1)
id.ys <- c(2,2,2,2,2)
ids.ys <- c(3,3,3,3,3)
iro <- c(1,3,2,1,"blue")
senshu <- c(1,1,1,1,1)
#  tenshu <- c(1,1)
tenshu <- c(19,3,3,2,22)
hutosa <- c(0.8,0.8,0.8,0.8,0.8)
is.leg <- 0

sen <- c(1,0,0,0,0)

phsiflag="phi"
xy <- "x"

name[1] <- paste(dirbase3,"/",pname3,"/pmf/1D/pmf_1D_UmbMD_vac_",TLns,"ns_",xy,"=",3,".txt",sep="")

name[2] <- paste(dirbase2,"/",pname2,"/nEX=",numEX2,"/freq=",TLbase2,"ps","/pmf/1D/pmf_py1",phsiflag,"@",3,sep="")

name[3] <- paste(dirbase1,"/tau=",tau[1],"/mZ=",mZ[1],"/TZ=",TZs[1],"/",pname1,"/nEX=",numEX1,"/fq=",TLbase1,"/pmf/1D/","pmf_pymbar_1_TAA=300_TCG=300_TZ=1000_KZAAo=1000_KZCGo=0_AA",phsiflag,"@",3,sep="")

name[4] <- paste(dirbase4,"/tau=",tauTAMD,"/TB=",TBTAMD,"/KZ=",KZTAMD,"/mZ=",mZTAMD,"/pmf_1D/pmf_1D_T=",TeTAMD,"_x=",3,".txt",sep="")

name[5] <- paste(dirbaseREUS,"/",pnameREUS,"/nEX_",numEXREUS,"/freq_",TLbaseREUS,"ps","/pmf/1D/pmf_1D_UmbMD_vac_nbx=20_nby=20_",xy,"=3",".txt",sep="")

label.x<-expression(paste(psi,"(radian,","/",pi,")"))

plmGeneralwsrange(data.names=name,
                  sd.names=name,
                  id.ys=id.ys,
                  ids.ys=ids.ys,
                  is.sen=sen,
                  label.size=0.5,axis.size=1.0,
                  iro=iro,axis.ft="F",is.header="T",
                  sdiro=iro,
                  xrange=xrange,yrange=yrange,
                  sdyrange=yrange,
                  warrow="T")
box(lwd=1.0)

txt <- "(a)"
text(0.75,10.5,txt,cex=1.4)
  
axis(2,yaxp=yrange.axis,lwd=1.0,cex.axis=1.0)


phsiflag="phi"
xy <- "x"

name[1] <- paste(dirbase3,"/",pname3,"/pmf/1D/pmf_1D_UmbMD_vac_",TLns,"ns_",xy,"=",6,".txt",sep="")

name[2] <- paste(dirbase2,"/",pname2,"/nEX=",numEX2,"/freq=",TLbase2,"ps","/pmf/1D/pmf_py1",phsiflag,"@",6,sep="")

name[3] <- paste(dirbase1,"/tau=",tau[1],"/mZ=",mZ[1],"/TZ=",TZs[1],"/",pname1,"/nEX=",numEX1,"/fq=",TLbase1,"/pmf/1D/","pmf_pymbar_1_TAA=300_TCG=300_TZ=1000_KZAAo=1000_KZCGo=0_AA",phsiflag,"@",6,sep="")

name[4] <- paste(dirbase4,"/tau=",tauTAMD,"/TB=",TBTAMD,"/KZ=",KZTAMD,"/mZ=",mZTAMD,"/pmf_1D/pmf_1D_T=",TeTAMD,"_x=",6,".txt",sep="")

name[5] <- paste(dirbaseREUS,"/",pnameREUS,"/nEX_",numEXREUS,"/freq_",TLbaseREUS,"ps","/pmf/1D/pmf_1D_UmbMD_vac_nbx=20_nby=20_",xy,"=6",".txt",sep="")

#cat(name)

fact.x <- 1/pi
ave.x <- 1

label.x<-expression(paste(psi,"(radian,","/",pi,")"))

xrange <- c(-1.0,1.0,4)
xrange.axis <- c(-1.0,1.0,4)
yrange <- c(0,12,6)
yrange.axis <- c(0,10,5)

id.ys <- c(2,2,2,2,2)
ids.ys <- c(3,3,3,3,3)

plmGeneralwsrange(data.names=name,
                  sd.names=name,
                  id.ys=id.ys,
                  ids.ys=ids.ys,
                  is.sen=sen,
                  label.size=0.5,axis.size=1.0,
                  iro=iro,axis.ft="F",is.header="T",
                  sdiro=iro,
                  xrange=xrange,yrange=yrange,
                  sdyrange=yrange,
                  warrow="T")
box(lwd=1.0)

txt <- "(b)"
text(0.75,10.5,txt,cex=1.4)

axis(2,yaxp=yrange.axis,lwd=1.0,cex.axis=1.0)
mtext(outer=T,label.y,side=2,line=2.0,cex=1.0)

axis(1,xaxp=xrange.axis,lwd=1.0,cex.axis=1.0)
mtext(label.x,side=1,line=2.0,cex=1.0)

dev.off()
EOF

Rscript ${paraRfil}; echo ${paraRfil}
