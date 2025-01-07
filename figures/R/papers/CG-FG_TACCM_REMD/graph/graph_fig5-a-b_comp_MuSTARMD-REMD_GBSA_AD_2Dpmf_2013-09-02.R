TAA_MuSTARMD<-"300"
TCG_MuSTARMD<-"300"

TZ_MuSTARMD <- "1000"

tau_MuSTARMD <- "1.0"

mZ_MuSTARMD <-"100.00"  

pname_MuSTARMD <- "SB_KZMAX=1000_NR8_woeljd0.001_2013-08-28"

nEX_MuSTARMD<-"1000"

fq_MuSTARMD<-"10"

TLbase_MuSTARMD<-"10"

width_MuSTARMD <- "0.3"

WVflag_MuSTARMD <- "InW"
pname_TREMD <-"f300t500_rep4"

numEX_TREMD <-"100000"

TLbase_TREMD <-"1"

ff_TREMD <-"ff99SB"

fact.x=1
fact.y=1

level <- seq(0,30,2)

name.title <- NULL

MyColor <- function(n,alpha=1)
{			
    if ((n <- as.integer(n[1L])) > 0) {
        j <- n%/%3		
        k <- n%/%3		
        i <- n - j - k		
        c(if (i > 0) hsv(h = seq.int(from = 40/60, to = 25/60,
            length.out = i), alpha = alpha), if (j > 0) hsv(h = seq.int(from = 23/60,
            to = 11/60, length.out = j), alpha = alpha), if (k > 
            0) hsv(h = seq.int(from = 8/60, to = 0/60, length.out = k-1),
            alpha = alpha, s = seq.int(from = 1, to = 0.9, length.out = k-1),
            v = 1),hsv(0,0,1))
    }			
    else character(0L)	
}			

file.name <- "~/papers/CG-FG_TACCM_REMD/eps/fig5-a-b_comp_MuSTARMD-TREMD_GBSA_AD_2Dpmf_2013-09-04.eps"
#tiff(file.name,width=720,height=300)
#postscript(file.name,width=3.6,height=4.3,horizontal=FALSE,onefile=FALSE,paper="special")
postscript(file.name,width=3.6,height=4.5,horizontal=FALSE,onefile=FALSE,paper="special")

xrange <- c(-1.0,1.0,4)
xrange.axis <- c(-1.0,1.0,4)
yrange <- c(-1.0,1.0,4)
yrange.axis <- c(-1.0,1.0,4)

par(mar = c(0.0,0.0,0.0,0.0) )
par(oma = c(0,0,0,0) )
source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange_worwoaxis.R")
source("~/papers/CG-FG_TACCM_REMD/FillConBar.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap_wrange_worwoaxis.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConBar_wmai.R")

nf <- layout(matrix(c(1,3,2,3),2,2,byrow=TRUE),c(3,1),c(20,23))

par(cex.axis=1.0)
par(cex.lab=1.0)

label.x <- expression(paste(phi,"(radian/",pi,")"))
#label.y <- expression(paste(,psi,"(radian/",pi,")"))

title=name.title
	
dir0 <- "/home/yamamori/calspa/InV2InW/AD"
dirbase <- paste(dir0,"/e_GBSA_2012-08-06/e_CG-FG_NH_2012-07-27/",sep="")
AACG <- "AA"
KZAAo <- "1000"
KZCGo <- "0"

name <- paste(dirbase,"/tau=",tau_MuSTARMD,"/mZ=",mZ_MuSTARMD,"/TZ=",TZ_MuSTARMD,"/",pname_MuSTARMD,"/nEX=",nEX_MuSTARMD,"/fq=",fq_MuSTARMD,"/pmf/pmf_TAA=300_TCG_300_TZ_1000_0.3_0_1000_CG_InW",sep="")
cat(name)

label.y <- expression(paste("ian/",pi,")","                                 "))

felwFillConMapwrangeworwoaxis(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,mai=c(0,0.8,0.5,0.1),xaflag="f",norm="T")

text(0.75,0.8,"(a)",cex=1.4)

dir0 <- "~/calspa/refcalc/REMD/AD"
dirbase <- paste(dir0,"/s_REGB_2012-07-23_",ff_TREMD,sep="")
name <- paste(dirbase,"/",pname_TREMD,"/nEX=",numEX_TREMD,"/freq=",TLbase_TREMD,"ps","/pmf/pmf_ADg",sep="")
cat(name)

label.y <- expression(paste("                                ",psi,"(rad"))

felwFillConMapwrangeworwoaxis(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,mai=c(0.8,0.8,0,0.1),norm="T")

text(0.75,0.8,"(b)",cex=1.4)

felwFillConBarwmai(name,label.x=label.x,label.y=label.y,level=level,mai=c(0.8,0.5,0.5,0.2))

