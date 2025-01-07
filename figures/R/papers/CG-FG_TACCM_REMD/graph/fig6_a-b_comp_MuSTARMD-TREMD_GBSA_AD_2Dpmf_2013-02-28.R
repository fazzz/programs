TAA_MuSTARMD<-"300"
TCG_MuSTARMD<-"300"

TZ_MuSTARMD <- "750"

tau_MuSTARMD <- "1.0"

mZ_MuSTARMD <-"50000.00"  

pname_MuSTARMD <- "T1_wrefd0.1"

nEX_MuSTARMD<-"1000"

fq_MuSTARMD<-"10"

TLbase_MuSTARMD<-"10"

width_MuSTARMD <- "0.3"

WVflag_MuSTARMD <- "InW"
pname_TREMD <-"f300t400"

numEX_TREMD <-"100000"

TLbase_TREMD <-"1"

ff_TREMD <-"ff99"

level <- seq(0,20,2)

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

file.name <- "~/papers/CG-FG_TACCM_REMD/eps/fig6_a-b_comp_MuSTARMD-TREMD_GBSA_AD_2Dpmf_2013-02-28.eps"
#tiff(file.name,width=720,height=300)
postscript(file.name,width=3.5,height=2.6,horizontal=FALSE,onefile=FALSE,paper="special")

xrange <- c(-3,3,6)
yrange <- c(-3,3,6)

par(oma = c(0,0,0,0) )
source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange.R")
source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange_wy.R")
source("~/papers/CG-FG_TACCM_REMD/FillConBar.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap_wrange.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap_wrange_wy.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConBar.R")

nf <- layout(matrix(c(1,2,3),1,3,byrow=TRUE),c(3.2,3,1),c(1))

par(mar =c(5,5,2,2) )
par(cex.axis=2.0)
par(cex.lab=2.0)
label.x <- expression(paste(phi))
label.y <- expression(paste(psi))

title=name.title
	
dir0 <- "/home/yamamori/calspa/InV2InW/AD"
dirbase <- paste(dir0,"/e_GBSA_2012-08-06/e_CG-FG_NH_2012-07-27/",sep="")
AACG <- "AA"
KZAAo <- "1000"
KZCGo <- "0"

name <- paste(dirbase,"/tau=",tau_MuSTARMD,"/mZ=",mZ_MuSTARMD,"/TZ=",TZ_MuSTARMD,"/",pname_MuSTARMD,"/nEX=",nEX_MuSTARMD,"/fq=",fq_MuSTARMD,"/pmf/pmf_TAA=",TAA_MuSTARMD,"_TCG_",TCG_MuSTARMD,"_TZ_",TZ_MuSTARMD,"_",width_MuSTARMD,"_",KZAAo,"_",KZCGo,"_",AACG,"_",WVflag_MuSTARMD,sep="")
cat(name)
felwFillConMapwrangewy(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,mar=c(5,5,1,0))

dir0 <- "~/calspa/refcalc/REMD/AD"
dirbase <- paste(dir0,"/s_REGB_2012-07-23_",ff_TREMD,sep="")
name <- paste(dirbase,"/",pname_TREMD,"/nEX=",numEX_TREMD,"/freq=",TLbase_TREMD,"ps","/pmf/pmf_ADg_radi",sep="")
cat(name)
felwFillConMapwrange(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,mar=c(5,0,1,0))

felwFillConBar(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

