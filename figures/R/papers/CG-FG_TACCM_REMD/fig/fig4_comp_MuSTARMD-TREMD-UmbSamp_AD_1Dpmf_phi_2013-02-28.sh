#!~/bin/sh

phsiflag=phi

parafile1=~/calspa/TACCM_CGAAREMD/AD/para/2012-11-02/para_e_CG-FG_NR=4_TZ=700_fq=10ps_99SB_KZmax=5000_weljd0.001_mZ=100.sh

parafile2=~/calspa/refcalc/REMD/AD/para/para_s_REMD_vac_ff99SB_2012-07-24.sh
parafileUmbSam=~/calspa/refcalc/UmbSam/AD/para/para_s_vac_ff99SB_K=10_2012-11-12.sh
height=20 
width1=1 
numuene=4 

dir=~/calspa/TACCM_CGAAREMD/AD
dirOUT=~/papers/CG-FG_TACCM_REMD

source ${parafile1}

height=${height} 

TZ=${TZs[1]}

AACGflag=CG

paraRfil=${dirOUT}/graph/para_comp_MuSTARMD-TREMD-UmbSamp_AD_1Dpmf_2013-01-14.R

cat <<EOF > ${paraRfil}
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

phsiflag <- "${phsiflag}"
phsi <- ${phsi}
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

name.out <- "${dirOUT}/tiff/fig4_comp_MuSTARMD-TREMD-UmbSamp_AD_1Dpmf_phi_2013-02-28"

source("~/Rspa/plmGeneral_wsrange.R")

dir01 <- "/home/yamamori/calspa/TACCM_CGAAREMD/AD"
dirbase1 <- paste(dir01,"/e_CG-FG_NH_2012-07-27",sep="")

dir02 <- "~/calspa/refcalc/REMD/AD"
dirbase2 <- paste(dir02,"/s_REVAC_2012-07-23_",ff,sep="")

dir03 <- "~/calspa/refcalc/UmbSam/AD"
dirbase3 <- paste(dir03,"/s_UmbSam_vac_2012-11-12_",ff2,sep="")

nTZs<-1

title=name.title

if ( phsiflag == "psi" )
  label.x<-expression(paste(phi))
if ( phsiflag == "phi" )
  label.x<-expression(paste(psi))

label.y="pmf"

xrange <- c(-3.14,3.14,4)
xrange.axis <- c(-2,3,5)
yrange.axis <- c(0,height,10)
height.ext <- height * 1.2
yrange <- c(0,height.ext,10)

fact.x <- 1
fact.y <- 1
ave.x <- 1

file.name <- paste(name.out,".tiff",sep='')
cat(file.name,'\n')
#tiff(file.name,width=800,height=700)
tiff(file.name,width=600,height=525)

par(mfrow=c(2,1))
par(mar=c(0.0,0.0,0.0,0.0))
par(oma=c(7.5,6.0,2.0,2.0))

s <- 1
name <- NULL

id.xs <- c(1,1,1)
id.ys <- c(2,2,2)
ids.ys <- c(3,3,3)
iro <- c(2,3,1)
senshu <- c(1,1,1)
#  tenshu <- c(1,1)
tenshu <- c(3,3,19)
hutosa <- c(1,1,1)
is.leg <- 0

if ( phsiflag=="phi" ) {
  xy <- "x"
}
if ( phsiflag=="psi" ) {
  xy <- "y"
}

name[1] <- paste(dirbase1,"/tau=",tau[1],"/mZ=",mZ[1],"/TZ=",TZs[1],"/",pname1,"/nEX=",numEX1,"/fq=",TLbase1,"/pmf/1D/pmf_pymbar_1_TAA=",TAA,"_TCG=",TCG,"_TZ=",TZs[1],"_KZAAo=",KZAAo,"_KZCGo=",KZCGo,"_",AACG,phsiflag,"@",6,sep="")
name[2] <- paste(dirbase2,"/",pname2,"/nEX=",numEX2,"/freq=",TLbase2,"ps","/pmf/1D/pmf_p1",phsiflag,"@",6,sep="")
name[3] <- paste(dirbase3,"/",pname3,"/pmf/1D/pmf_1D_UmbMD_vac_",TLns,"ns_",xy,"=",6,".txt",sep="")

plmGeneralwsrange(data.names=name,
                  sd.names=name,
                  id.ys=id.ys,
                  ids.ys=ids.ys,
                  is.sen=rep(0,20),
                  label.size=0.5,axis.size=2.0,
                  iro=iro,axis.ft="F",is.header="T",
                  sdiro=iro,
                  xrange=xrange,yrange=yrange,
                  sdyrange=yrange,
                  warrow="T")
box(lwd=2.0)
  
axis(2,yaxp=yrange.axis,lwd=2.0,cex.axis=1.5)
mtext(outer=T,label.y,side=2,line=4.0,cex=1.5)

name[1] <- paste(dirbase1,"/tau=",tau[1],"/mZ=",mZ[1],"/TZ=",TZs[1],"/",pname1,"/nEX=",numEX1,"/fq=",TLbase1,"/pmf/1D/pmf_pymbar_1_TAA=",TAA,"_TCG=",TCG,"_TZ=",TZs[1],"_KZAAo=",KZAAo,"_KZCGo=",KZCGo,"_",AACG,phsiflag,"@",12,sep="")
name[2] <- paste(dirbase2,"/",pname2,"/nEX=",numEX2,"/freq=",TLbase2,"ps","/pmf/1D/pmf_p1",phsiflag,"@",12,sep="")
name[3] <- paste(dirbase3,"/",pname3,"/pmf/1D/pmf_1D_UmbMD_vac_",TLns,"ns_",xy,"=",12,".txt",sep="")

cat(name)

plmGeneralwsrange(data.names=name,
                  sd.names=name,
                  id.ys=id.ys,
                  ids.ys=ids.ys,
                  is.sen=rep(0,20),
                  label.size=0.5,axis.size=2.0,
                  iro=iro,axis.ft="F",is.header="T",
                  sdiro=iro,
                  xrange=xrange,yrange=yrange,
                  sdyrange=yrange,
                  warrow="T")
box(lwd=2.0)
  
axis(2,yaxp=yrange.axis,lwd=2.0,cex.axis=1.5)
mtext(outer=T,label.y,side=2,line=4.0,cex=1.5)
  
axis(1,xaxp=xrange.axis,lwd=2.0,cex.axis=1.5)
mtext(outer=T,label.x,side=1,line=4.0,cex=1.5)

dev.off()
EOF

Rscript ${paraRfil}; echo ${paraRfil}



