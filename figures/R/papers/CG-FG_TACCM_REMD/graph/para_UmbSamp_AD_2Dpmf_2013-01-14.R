fact.x <- 1
fact.y <- 1
fact.p <- 1.0/pi

dir    <- "/home/yamamori/calspa/refcalc/UmbSam/AD"
dirout <- "/home/yamamori/papers/CG-FG_TACCM_REMD"

title <- NULL

#label.x <- expression(paste(phi,"(radian/",pi,")"))
#label.y <- expression(paste(psi,"(radian/",pi,")"))

label.x <- expression(paste(""))
label.y <- expression(paste(""))

name.out <- paste(dirout,"/eps/","fig3-a-g_comp_MuSTARMD-REMD-REUS-TAMD-CMD_AD_2Dpmf_2013-08-20",sep='')
   
#level <- seq(0.0,10.0,1.0)
level <- seq(0.0,30.0,2.0)

pdatanamedummy <- paste(dirout,"/fig/points_dummy.txt",sep='')
pdataname <- paste(dirout,"/fig/points.txt",sep='')
ldataname <- paste(dirout,"/fig/lines.txt",sep='')
ldataname2 <- paste(dirout,"/fig/lines2.txt",sep='')
ldataname3 <- paste(dirout,"/fig/lines3.txt",sep='')

xrange <- c(-1.0,1.0,4)
xrange.axis <- c(-1.0,1.0,4)
yrange <- c(-1.0,1.0,4)
yrange.axis <- c(-0.5,1.0,3)

file.name <- paste(name.out,'.eps',sep='')
postscript(file.name,width=5.3,height=6.0,horizontal=FALSE,onefile=FALSE,paper="special")

par(oma = c(0,0,0,0) )
source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange_worwoaxis.R")
source("~/papers/CG-FG_TACCM_REMD/FillConBar.R")

source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap_wrange_worwoaxis_wpoints.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConBar_wmai.R")

source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMapbox_wpoints_2.R")

nf <- layout(matrix(c(1,2,7,3,4,7,5,6,7),3,3,byrow=TRUE),c(22.5,16.5,8),c(20,15,23))

#par(cex.axis=1.0)
#par(cex.lab=1.0)

par(cex.axis=1.44)
par(cex.lab=1.44)

name <- paste("/home/yamamori/calspa/refcalc/UmbSam/AD/s_UmbSam_vac_2012-11-12_ff99SB/Umb_Nbin=12x12_K=10/pmf/pmf_UmbMD_vac_10ns.txt_2",sep='')

cat(name)

#label.x <- expression(paste(phi,"(radian/",pi,")"))
label.y <- ""

felwFillConMapwrangeworwoaxiswpoints(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,mai=c(0.0,0.8,0.5,0.0),xaflag="f",norm="T",pdata=pdatanamedummy,iro=c("white","white","white"))

text(0.75,0.8,"(a)",cex=1.6)

TAA_MuSTARMD<-"300"
TCG_MuSTARMD<-"300"

TZ_MuSTARMD <- "700"

tau_MuSTARMD <- "1.0"

mZ_MuSTARMD <-"100.00"  

pname_MuSTARMD <- "SB_KZMAX=500_NR4_woeljd0.001"

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
name <- paste(dirbase,"/tau=",tau_MuSTARMD,"/mZ=",mZ_MuSTARMD,"/TZ=",TZ_MuSTARMD,"/",pname_MuSTARMD,"/nEX=",numEX_MuSTARMD,"/fq=",fq_MuSTARMD,"/pmf/pmf_pymbar_TAA=300_TCG=300_TZ=700_KZAAo=500_KZCGo=0_AA_2013-05-26_20_2",sep="") 

xrange <- c(-1.0,1.0,4)
xrange.axis <- c(-1.0,1.0,4)
yrange <- c(-1.0,1.0,4)
yrange.axis <- c(-1.0,0.5,3)

cat(name)

#label.x <- expression(paste(phi,"(radian/",pi,")"))
label.y <- ""

felwFillConMapwrangeworwoaxiswpoints(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,mai=c(0.0,0.0,0.5,0.1),xaflag="f",yaflag="f",norm="T",pdata=pdatanamedummy,iro=c("white","white","white"))

text(0.75,0.8,"(b)",cex=1.6)

pname_TREMD <-"f300t600_rep4"

numEX_TREMD <-"100000"

TLbase_TREMD <-"1"

ff_TREMD <-"ff99SB"

dir0 <- "~/calspa/refcalc/REMD/AD"
dirbase <- paste(dir0,"/s_REVAC_2012-07-23_",ff_TREMD,sep="")
name <- paste(dirbase,"/",pname_TREMD,"/nEX=",numEX_TREMD,"/freq=",TLbase_TREMD,"ps","/pmf/pmf_py3",sep="")

xrange <- c(-1.0,1.0,4)
xrange.axis <- c(-1.0,1.0,4)
yrange <- c(-1.0,1.0,4)
yrange.axis <- c(-1.0,1.0,4)

#label.x <- expression(paste(phi,"(radian/",pi,")"))
#label.y <- expression(paste(psi,"(radian/",pi,")"))

cat(name)

felwFillConMapwrangeworwoaxiswpoints(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,mai=c(0.0,0.8,0.0,0.0),xaflag="f",norm="T",pdata=pdatanamedummy,iro=c("white","white","white"))

text(0.75,0.8,"(c)",cex=1.6)

pname_REUS <-"UmbAt_Nbin_4x4_K_1"

numEX_REUS <-"250"

TLbase_REUS <-"10"

ff_REUS <-"ff99SB"

dir0 <- "~/calspa/refcalc/REUS/AD"
dirbase <- paste(dir0,"/s_REUSVAC_2013-08-12_",ff_REUS,sep="")
name <- paste(dirbase,"/",pname_REUS,"/nEX_",numEX_REUS,"/freq_",TLbase_REUS,"ps","/pmf/pmf_UmbMD_vac_nbx=40_nby=40.txt_2",sep="")

xrange <- c(-1.0,1.0,4)
xrange.axis <- c(-1.0,1.0,4)
yrange <- c(-1.0,1.0,4)
yrange.axis <- c(-1.0,1.0,4)

cat(name)

felwFillConMapwrangeworwoaxiswpoints(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,mai=c(0.0,0.0,0.0,0.1),xaflag="f",yaflag="f",norm="T",pdata=pdatanamedummy,iro=c("white","white","white"))

text(0.75,0.8,"(d)",cex=1.6)

T_TAMD<-"300"

TB_TAMD<-"750"

tau_TAMD<-"1.0"

KZ_TAMD<-"1000"

mZ_TAMD<-"100.00"

width_TAMD <- "0.3"

ff_TAMD <- "99SB"

dir0 <- "/home/yamamori/calspa/TACCM/AD"
dirbase <- paste(dir0,"/e_TACCM_NH_2012-07-31_",ff_TAMD,sep="")
name <- paste(dirbase,"/tau=",tau_TAMD,"/TB=",TB_TAMD,"/KZ=",KZ_TAMD,"/mZ=",mZ_TAMD,"/pmf/pmf_T=",T_TAMD,".Zhist_0.2",sep="")

xrange <- c(-1.0,0.5,3)
xrange.axis <- c(-1.0,0.5,3)
yrange <- c(-1.0,0.5,3)
yrange.axis <- c(-1.0,0.5,3)

#label.x <- expression(paste(phi,"(radian/",pi,")"))
#label.y <- ""

cat(name)

felwFillConMapwrangeworwoaxiswpoints(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,mai=c(0.8,0.8,0.0,0.0),norm="T",pdata=pdatanamedummy,iro=c("white","white","black"))

text(0.75,0.8,"(e)",cex=1.6)

pname_CMD<-"UmbAt_Nbin_4x4_K_1"

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

#label.x <- expression(paste("                                 ",phi,"(rad"))
#label.y <- ""

cat(name)

felwFillConMapwrangeworwoaxisboxwpoints(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,mai=c(0.8,0.0,0.0,0.1),norm="T",yaflag="f",pdata=pdataname,ldata=ldataname,ldata2=ldataname2,ldata3=ldataname3,iro=c("pink","pink","pink","red","red","red"))

text(-0.6,0.8,"C5",cex=1.2,col="white")
text(-0.5,0.20,"C7eq",cex=1.2,col="white")
text(0.6,-0.2,"C7ax",cex=1.2)
text(0.75,0.8,"(f)",cex=1.6)


felwFillConBarwmai(name,label.x=label.x,label.y=label.y,level=level,mai=c(4.2,0.5,0.5,0.2))

