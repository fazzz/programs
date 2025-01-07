TAA_MuSTARMD<-"300"
TCG_MuSTARMD<-"300"

TZ_MuSTARMD <- "750"

tau_MuSTARMD <- "1.0"

mZ_MuSTARMD <-"50000.00"  

pname_MuSTARMD <- "99SB_T1_wrefd0.1"

numEX_MuSTARMD<-"1000"

fq_MuSTARMD<-"10"

TLbase_MuSTARMD<-"10"

width_MuSTARMD <- "0.3"

pname_TREMD <-"f300t400"

numEX_TREMD <-"100000"

TLbase_TREMD <-"1"

ff_TREMD <-"ff99SB"

T_TAMD<-"300"

TB_TAMD<-"750"

tau_TAMD<-"1.0"

KZ_TAMD<-"1000"

mZ_TAMD<-"50000.00"

width_TAMD <- "0.3"

ff_TAMD <- "99SB"

pname_CMD<-"f300t400"

T_CMD<-"300"

ff_CMD<-"ff99SB"

width_CMD <- "0.2"

dir    <- "/home/yamamori/calspa/refcalc/UmbSam/AD"
dirout <- "/home/yamamori/papers/CG-FG_TACCM_REMD"

title <- NULL

label.x <- expression(paste(phi,"(radian/",pi,")"))
label.y <- expression(paste(psi,"(radian/",pi,")"))

name.out <- paste(dirout,"/eps/","fig3_a-b-c-d-e-f_comp_MuSTARMD-other_methods_AD_2Dpmf_2013-05-08",sep='')

level <- seq(0,20,2)

xrange <- c(-1.0,1.0,4)
yrange <- c(-1.0,1.0,4)
   
file.name <- paste(name.out,'.eps',sep='')
#tiff(file.name,width=,height=)
#postscript(file.name,width=5.6,height=5.8,horizontal=FALSE,onefile=FALSE,paper="special")
postscript(file.name,width=10,height=10,horizontal=FALSE,onefile=FALSE,paper="special")

par(oma = c(0,0,0,0) )
source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange_worwoaxis.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap_wrange_worwoaxis.R")
source("~/papers/CG-FG_TACCM_REMD/FillConBar.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConBar_wmai.R")

nf <- layout(matrix(c(1,2,7,3,4,7,5,6,7),3,3,byrow=TRUE),c(23,16,12),c(20,15,23))

fact.x <- 1
fact.y <- 1

par(cex.axis=0.8)
par(cex.lab=1.0)

name.pmf <- paste("/home/yamamori/calspa/refcalc/UmbSam/AD/s_UmbSam_vac_2012-11-12_ff99SB/Umb_Nbin=12x12_K=10/pmf/pmf_UmbMD_vac_10ns.txt",sep='')

cat(name.pmf,"\n")

#felwFillConMapwrangeworwoaxis(name.pmf,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,mai=c(0.0,0.8,0.5,0.0),xaflag="f",norm="T")

felwFillConMapwrangeworwoaxis(name.pmf,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xaflag="f",norm="T")

dir0 <- "/home/yamamori/calspa/TACCM_CGAAREMD/AD"
dirbase <- paste(dir0,"/e_CG-FG_NH_2012-07-27",sep="")
AACG <- "CG"
KZAAo <- "0"
KZCGo <- "1000"
name <- paste(dirbase,"/tau=",tau_MuSTARMD,"/mZ=",mZ_MuSTARMD,"/TZ=",TZ_MuSTARMD,"/",pname_MuSTARMD,"/nEX=",numEX_MuSTARMD,"/fq=",fq_MuSTARMD,"/pmf/pmf_TAA=",TAA_MuSTARMD,"_TCG_",TCG_MuSTARMD,"_TZ_",TZ_MuSTARMD,"_",width_MuSTARMD,"_",KZAAo,"_",KZCGo,"_",AACG,"_4_2012-08-21",sep="") 

#felwFillConMapwrangeworwoaxis(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,mai=c(0.0,0.0,0.5,0.1),xaflag="f",yaflag="f",,norm="T")

felwFillConMapwrangeworwoaxis(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xaflag="f",yaflag="f",,norm="T")

dir0 <- "~/calspa/refcalc/REMD/AD"
dirbase <- paste(dir0,"/s_REVAC_2012-07-23_",ff_TREMD,sep="")
name <- paste(dirbase,"/",pname_TREMD,"/nEX=",numEX_TREMD,"/freq=",TLbase_TREMD,"ps","/pmf/pmf_ADv",sep="")

#felwFillConMapwrangeworwoaxis(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,mai=c(0.0,0.8,0.0,0.0),xaflag="f",norm="T")

felwFillConMapwrangeworwoaxis(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xaflag="f",norm="T")

dir0 <- "/home/yamamori/calspa/TACCM/AD"
dirbase <- paste(dir0,"/e_TACCM_NH_2012-07-31_",ff_TAMD,sep="")
name <- paste(dirbase,"/tau=",tau_TAMD,"/TB=",TB_TAMD,"/KZ=",KZ_TAMD,"/mZ=",mZ_TAMD,"/pmf/pmf_T=",T_TAMD,".Zhist",sep="")

#felwFillConMapwrangeworwoaxis(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,mai=c(0.0,0.0,0.0,0.1),xaflag="f",yaflag="f",norm="T")

felwFillConMapwrangeworwoaxis(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xaflag="f",yaflag="f",norm="T")

dir0 <- "~/calspa/refcalc/CMD/AD"
dirbase <- paste(dir0,"/s_CVAC_2012-08-08_",ff_CMD,sep="")
name <- paste(dirbase,"/anl/pmf_ADv_T",T_CMD,"_",width_CMD,sep="")

#felwFillConMapwrangeworwoaxis(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,mai=c(0.8,0.8,0.0,0.0),xaflag="f",norm="T")

felwFillConMapwrangeworwoaxis(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xaflag="f",norm="T")

#felwFillConBarwmai(name.pmf,label.x=label.x,label.y=label.y,level=level,mai=c(0.8,0.5,0.5,0.2))

#dev.off()
