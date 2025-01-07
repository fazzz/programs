fact.x <- 1
fact.y <- 1
fact.p <- 1.0/pi

dir    <- "/home/yamamori/calspa/refcalc/UmbSam/AD"
dirout <- "/home/yamamori/thesis/MuSTARMD"

title <- NULL

label.x <- expression(paste(phi,"(radian/",pi,")"))
label.y <- expression(paste(psi,"(radian/",pi,")"))

name.out <- paste(dirout,"/eps/","fig_Chapt6-2-a-b-c-d-e-f_comp_MuSTARMD-other_methods_AD_2Dpmf",sep='')
   
#level <- seq(0.0,10.0,1.0)
level <- seq(0.0,20.0,2.0)

pdatanamedummy <- paste(dirout,"/fig/points_dummy.txt",sep='')
pdataname <- paste(dirout,"/fig/points.txt",sep='')
ldataname <- paste(dirout,"/fig/lines.txt",sep='')
ldataname2 <- paste(dirout,"/fig/lines2.txt",sep='')
ldataname3 <- paste(dirout,"/fig/lines3.txt",sep='')

xrange <- c(-1.0,1.0,4)
xrange.axis <- c(-1.0,1.0,4)
yrange <- c(-1.0,1.0,4)
yrange.axis <- c(-1.0,1.0,4)

file.name <- paste(name.out,'.eps',sep='')
postscript(file.name,width=5.3,height=6.0,horizontal=FALSE,onefile=FALSE,paper="special")

par(oma = c(0,0,0,0) )
source("~/thesis/MuSTARMD/FillConMap_wrange_worwoaxis.R")
source("~/thesis/MuSTARMD/FillConBar.R")

source("~/thesis/MuSTARMD/fel_wFillConMap_wrange_worwoaxis_wpoints.R")
source("~/thesis/MuSTARMD/fel_wFillConBar_wmai.R")

source("~/thesis/MuSTARMD/fel_wFillConMapbox_wpoints.R")

nf <- layout(matrix(c(1,2,7,3,4,7,5,6,7),3,3,byrow=TRUE),c(23,16,8),c(20,15,23))

#par(cex.axis=1.0)
#par(cex.lab=1.0)

par(cex.axis=1.44)
par(cex.lab=1.44)

name <- paste("/home/yamamori/calspa/refcalc/UmbSam/AD/s_UmbSam_vac_2012-11-12_ff99SB/Umb_Nbin=12x12_K=10/pmf/pmf_UmbMD_vac_10ns.txt_2",sep='')

cat(name)

label.x <- expression(paste(phi,"(radian/",pi,")"))
label.y <- ""

felwFillConMapwrangeworwoaxiswpoints(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,mai=c(0.0,0.8,0.5,0.0),xaflag="f",norm="T",pdata=pdatanamedummy,iro=c("white","white","white"))

text(0.75,0.8,"(a)",cex=1.6)

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

xrange <- c(-1.0,1.0,4)
xrange.axis <- c(-1.0,1.0,4)
yrange <- c(-1.0,1.0,4)
yrange.axis <- c(-1.0,0.5,3)

cat(name)

label.x <- expression(paste(phi,"(radian/",pi,")"))
label.y <- ""

felwFillConMapwrangeworwoaxiswpoints(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,mai=c(0.0,0.0,0.5,0.1),xaflag="f",yaflag="f",norm="T",pdata=pdatanamedummy,iro=c("white","white","white"))

text(0.75,0.8,"(b)",cex=1.6)

pname_TREMD <-"f300t400"

numEX_TREMD <-"100000"

TLbase_TREMD <-"1"

ff_TREMD <-"ff99SB"

dir0 <- "~/calspa/refcalc/REMD/AD"
dirbase <- paste(dir0,"/s_REVAC_2012-07-23_",ff_TREMD,sep="")
name <- paste(dirbase,"/",pname_TREMD,"/nEX=",numEX_TREMD,"/freq=",TLbase_TREMD,"ps","/pmf/pmf_ADv",sep="")

xrange <- c(-1.0,1.0,4)
xrange.axis <- c(-1.0,1.0,4)
yrange <- c(-1.0,1.0,4)
yrange.axis <- c(-1.0,1.0,4)

label.x <- expression(paste(phi,"(radian/",pi,")"))
label.y <- expression(paste(psi,"(radian/",pi,")"))

felwFillConMapwrangeworwoaxiswpoints(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,mai=c(0.0,0.8,0.0,0.0),xaflag="f",norm="T",pdata=pdatanamedummy,iro=c("white","white","white"))

text(0.75,0.8,"(c)",cex=1.6)

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

xrange <- c(-1.0,1.0,4)
xrange.axis <- c(-1.0,1.0,4)
yrange <- c(-1.0,1.0,4)
yrange.axis <- c(-1.0,1.0,4)

label.x <- expression(paste(phi,"(radian/",pi,")"))
label.y <- ""

felwFillConMapwrangeworwoaxiswpoints(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,mai=c(0.0,0.0,0.0,0.1),xaflag="f",yaflag="f",norm="T",pdata=pdatanamedummy,iro=c("white","white","white"))

text(0.75,0.8,"(d)",cex=1.6)

pname_CMD<-"f300t400"

T_CMD<-"300"

ff_CMD<-"ff99SB"

width_CMD <- "0.2"

dir0 <- "~/calspa/refcalc/CMD/AD"
dirbase <- paste(dir0,"/s_CVAC_2012-08-08_",ff_CMD,sep="")
name <- paste(dirbase,"/anl/pmf_ADv_T",T_CMD,"_",width_CMD,sep="")

xrange <- c(-1.0,1.0,4)
xrange.axis <- c(-1.0,1.0,4)
yrange <- c(-1.0,1.0,4)
yrange.axis <- c(-1.0,0.5,3)

label.x <- expression(paste("                                 ",phi,"(rad"))
label.y <- ""

felwFillConMapwrangeworwoaxiswpoints(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,mai=c(0.8,0.8,0.0,0.0),norm="T",pdata=pdatanamedummy,iro=c("white","white","black"))

text(0.75,0.8,"(e)",cex=1.6)


xrange <- c(-1.0,1.0,4)
xrange.axis <- c(-0.5,1.0,3)
yrange <- c(-1.0,1.0,4)
yrange.axis <- c(-1.0,1.0,4)

label.x <- expression(paste("ian/",pi,")","                                "))
label.y <- ""

felwFillConMapwrangeworwoaxisboxwpoints(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,mai=c(0.8,0.0,0.0,0.1),yaflag="f",norm
="T",pdata=pdataname,ldata=ldataname,ldata2=ldataname2,ldata3=ldataname3,iro=c("black","black","black","red","red","red"))

text(-0.6,0.8,"C5",cex=1.2)
text(-0.2,0.45,"C7eq",cex=1.2)
text(0.6,-0.2,"C7ax",cex=1.2)

text(0.75,0.8,"(f)",cex=1.6)


felwFillConBarwmai(name,label.x=label.x,label.y=label.y,level=level,mai=c(0.8,0.5,0.5,0.2))

