fact.y=1

TAA<-"300"
TCG<-"300"

TZs <- c( "1000" )
nTZs <- length(TZs)

tau<-c(  "1.0" )
ntau<-length(tau)

mZ<-c(  "100.00" )
nmZ<-length(mZ)

nKZAA<-8
nKZCG<-8
numRE<-8

pname1<-"SB_KZMAX=1000_NR8_woeljd0.001_2013-08-28"

numEX1<-"1000"

TLbase1<-"10"

AACG <- "CG"
width <- "1"

KZAAo <- "0"
KZCGo <- "1000"

numuene <- "4"

height <- 20

phsiflag <- "phi"
phsi <- "phi"
name.title <- paste(AACG,sep="")


ff<-"ff99SB"

pname2<-"f300t1200"

numEX2<-"100000"

TLbase2<-"1"
TLbaseREUS<-"10"

numEXREUS<-"500"

numREUS<-"16"

pnameREUS<-"UmbAt_Nbin_4x4_K_0.5"
ff2<-"ff99SB"

pname3<-"Umb_Nbin=12x12_K=10"

TLns<-"10"

ffTAMD <-"99SB_80ns"
tauTAMD <-"1.0"
TBTAMD <- 800
KZTAMD <- "1000"
mZTAMD <- "50000.00"
TeTAMD=300


name.out <- "/home/yamamori/papers/CG-FG_TACCM_REMD/eps/fig4-c_comp_MuSTARMD-TREMD-REUS-TAMD-UmbSamp_AD_1Dpmf_2013-09-04"

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

xrange <- c(0,10,3)
xrange.axis <- c(0,10,3)

yrange <- c(0,18,3)     
yrange.axis <- c(0,18,3)

fact.x <- 1/pi
#fact.y <- 1/pi
ave.x <- 1

file.name <- paste(name.out,".eps",sep='')
#cat(file.name,'\n')
#tiff(file.name,width=800,height=700)
#tiff(file.name,width=600,height=525)
postscript(file.name,width=2.7,height=2.0,horizontal=FALSE,onefile=FALSE,paper="special")

label.x="C7eq                              C7ax"

par(mfrow=c(1,1))
par(mai=c(0.0,0.0,0.0,0.0))
par(oma = c(2.0,3.0,1.0,0.6) )

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

yrange <- c(0,12,6)
yrange.axis <- c(0,12,6)

phsiflag="psi"
xy <- "y"

sen <- c(1,0,0,0,0)

name[1] <- paste("~/papers/CG-FG_TACCM_REMD/FEL_C7eq-C7ax_UmbSamMD.txt")
name[2] <- paste("~/papers/CG-FG_TACCM_REMD/FEL_C7eq-C7ax_REMD_2013-09-04.txt")
name[3] <- paste("~/papers/CG-FG_TACCM_REMD/FEL_C7eq-C7ax_MuSTARMD_2013-09-04.txt")
name[4] <- paste("~/papers/CG-FG_TACCM_REMD/FEL_C7eq-C7ax_TAMD_2013-09-04.txt")
name[5] <- paste("~/papers/CG-FG_TACCM_REMD/FEL_C7eq-C7ax_REUS_2013-09-04.txt")

cat(name)

fact.x <- 1
fact.y <- 1
ave.x <- 1

id.ys <- c(4,4,4,4,4)
ids.ys <- c(5,5,5,5,5)

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

txt <- "(c)"
#text(9,10,txt,cex=1.4)

#labelname<-c("","","","","","","","","","","")
#axis(1,xaxp=xrange,lwd=1.0,cex.axis=1.0)
#mtext(outer=T,label.x,side=1,line=0.5,cex=1.0)

axis(2,yaxp=yrange.axis,lwd=1.0,cex.axis=1.0)
mtext(outer=T,label.y,side=2,line=2.0,cex=1.0)

dev.off()
