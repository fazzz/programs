fact.x <- 1
fact.y <- 1
fact.p <- 1.0/pi

dir    <- ""
dirout <- "/home/yamamori/thesis/abstract"

title <- NULL

label.x <- expression(paste(""))
label.y <- expression(paste(""))

name.out <- paste(dirout,"/eps/","fig1_FEL_comp_MuSTARMD-others_AD_2013-06-03",sep='')
   
level <- seq(0.0,20.0,2.0)

file.name <- paste(name.out,'.eps',sep='')
postscript(file.name,width=3.2,height=5.8,horizontal=FALSE,onefile=FALSE,paper="special")

par(oma = c(0,0,0,0))

source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange_worwoaxis.R")
source("~/papers/CG-FG_TACCM_REMD/FillConBar.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap_wrange_worwoaxis.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConBar_wmai.R")

nf <- layout(matrix(c(1,4,2,4,3,4),3,2,byrow=TRUE),c(3,1),c(20,15,23))

par(cex.axis=1.44)
par(cex.lab=1.44)

TAA_MuSTARMD<-"300"
TCG_MuSTARMD<-"300"

TZ_MuSTARMD <- "750"

tau_MuSTARMD <- "1.0"

mZ_MuSTARMD <-"100.00"  

pname_MuSTARMD <- "SB_KZMAX=5000_NR8-2_woeljd0.001"

numEX_MuSTARMD<-"1000"

fq_MuSTARMD<-"10"

TLbase_MuSTARMD<-"10"

width_MuSTARMD <- "0.3"

dir0 <- "/home/yamamori/calspa/TACCM_CGAAREMD/AD"
dirbase <- paste(dir0,"/e_CG-FG_NH_2012-07-27",sep="")
AACG <- "CG"
KZAAo <- "0"
#KZCGo <- "1000"
KZCGo <- "5000"
name <- paste(dirbase,"/tau=",tau_MuSTARMD,"/mZ=",mZ_MuSTARMD,"/TZ=",TZ_MuSTARMD,"/",pname_MuSTARMD,"/nEX=",numEX_MuSTARMD,"/fq=",fq_MuSTARMD,"/pmf/pmf_TAA=",TAA_MuSTARMD,"_TCG_",TCG_MuSTARMD,"_TZ_",TZ_MuSTARMD,"_",width_MuSTARMD,"_",KZAAo,"_",KZCGo,"_",AACG,"_4_2012-08-21",sep="") 
#name <- paste(dirbase,"/tau=",tau_MuSTARMD,"/mZ=",mZ_MuSTARMD,"/TZ=",TZ_MuSTARMD,"/",pname_MuSTARMD,"/nEX=",numEX_MuSTARMD,"/fq=",fq_MuSTARMD,"/pmf/pmf_pymbar_TAA=",TAA_MuSTARMD,"_TCG=",TCG_MuSTARMD,"_TZ=",TZ_MuSTARMD,"_KZAAo=",KZAAo,"_KZCGo=",KZCGo,"_",AACG,sep="") 

xrange <- c(-3.0,3.0,6)
xrange.axis <- c(-3.0,3.0,6)
yrange <- c(-3.0,3.0,6)
yrange.axis <- c(-2.0,3.0,5)

cat(name)

felwFillConMapwrangeworwoaxis(name,label.x=label.x,label.y=label.y,level=level,title=title,xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,mai=c(0.0,0.8,0.5,0.1),xaflag="f")

text(2.3,2.4,"(a)",cex=1.6)

pname_TREMD <-"f300t400"

numEX_TREMD <-"100000"

TLbase_TREMD <-"1"

ff_TREMD <-"ff99SB"

dir0 <- "~/calspa/refcalc/REMD/AD"
dirbase <- paste(dir0,"/s_REVAC_2012-07-23_",ff_TREMD,sep="")
name <- paste(dirbase,"/",pname_TREMD,"/nEX=",numEX_TREMD,"/freq=",TLbase_TREMD,"ps","/pmf/pmf_ADv",sep="")

xrange <- c(-3.0,3.0,5)
xrange.axis <- c(-3.0,3.0,5)
yrange <- c(-3.0,3.0,5)
yrange.axis <- c(-2.0,3.0,5)

#label.x <- expression(paste(phi,"(radian)"))
label.y <- expression(paste(psi,"(radian)"))

felwFillConMapwrangeworwoaxis(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,mai=c(0.0,0.8,0.0,0.1),xaflag="f")

text(2.3,2.4,"(b)",cex=1.6)

T_TAMD<-"300"

TB_TAMD<-"750"

tau_TAMD<-"1.0"

KZ_TAMD<-"1000"

mZ_TAMD<-"50000.00"

width_TAMD <- "0.3"

ff_TAMD <- "99SB"

dir0 <- "/home/yamamori/calspa/TACCM/AD"
dirbase <- paste(dir0,"/e_TACCM_NH_2012-07-31_",ff_TAMD,sep="")
name <- paste(dirbase,"/tau=",tau_TAMD,"/TB=",TB_TAMD,"/KZ=",KZ_TAMD,"/mZ=",mZ_TAMD,"/pmf/pmf_T=",T_TAMD,".Zhist",sep="")

xrange <- c(-3.0,3.0,6)
xrange.axis <- c(-3.0,3.0,6)
yrange <- c(-3.0,3.0,6)
yrange.axis <- c(-3.0,3.0,6)

label.x <- expression(paste(phi,"(radian)"))
label.y <- ""

felwFillConMapwrangeworwoaxis(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,mai=c(0.8,0.8,0.0,0.1))

text(2.3,2.4,"(c)",cex=1.6)


felwFillConBarwmai(name,label.x=label.x,label.y=label.y,level=level,mai=c(0.8,0.5,0.5,0.2))

