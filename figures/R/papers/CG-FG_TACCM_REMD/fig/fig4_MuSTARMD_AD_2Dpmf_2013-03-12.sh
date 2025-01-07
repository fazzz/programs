#!/bin/zsh

dir=~/papers/CG-FG_TACCM_REMD/

paraRfil=${dir}/graph/graph_MuSTARMD_AD_2Dpmf_2013-02-28.R

tiffwidth=800
tiffheight=600

parafile_MuSTARMD=~/calspa/TACCM_CGAAREMD/AD/para/para_e_CG-FG_NH_2012-07-27_wrefd0.1_TZ=750_fq=10ps_99SB.sh 
source ${parafile_MuSTARMD}
cat << eof > ${paraRfil}
TAA_MuSTARMD<-"${TAA}"
TCG_MuSTARMD<-"${TCG}"

TZ_MuSTARMD <- "${TZs[1]}"

tau_MuSTARMD <- "${tau[1]}"

mZ_MuSTARMD <-"${mZ[1]}"  

pname_MuSTARMD <- "${pname}"

numEX_MuSTARMD<-"${numEX}"

fq_MuSTARMD<-"${TLbase}"

TLbase_MuSTARMD<-"${TLbase}"

width_MuSTARMD <- "0.3"

eof

cat << eof >> ${paraRfil}
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

file.name <- "~/papers/CG-FG_TACCM_REMD/eps/fig4_MuSTARMD_AD_2Dpmf_2013-02-28.eps"
#tiff(file.name,width=${tiffwidth},height=${tiffheight})
postscript(file.name,width=3.0,height=2.25,horizontal=FALSE,onefile=FALSE,paper="special")

xrange <- c(-3,3,6)
yrange <- c(-3,3,6)

par(oma = c(0,0,0,0) )
source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange.R")
source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange_wy.R")
source("~/papers/CG-FG_TACCM_REMD/FillConBar.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap_wrange.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap_wrange_wy.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConBar.R")
						
nf <- layout(matrix(c(1,2),1,2,byrow=TRUE),c(3,1),c(1))

#par(mar =c(5,5,2,2) )
#par(cex.axis=2.0)
#par(cex.lab=2.0)
par(mar =c(2.5,2.5,1,1) )
par(cex.axis=1.0)
par(cex.lab=1.0)
label.x <- expression(paste(phi))
label.y <- expression(paste(psi))

title=name.title
	
dir0 <- "/home/yamamori/calspa/TACCM_CGAAREMD/AD"
dirbase <- paste(dir0,"/e_CG-FG_NH_2012-07-27",sep="")
AACG <- "CG"
KZAAo <- "0"
KZCGo <- "1000"
name <- paste(dirbase,"/tau=",tau_MuSTARMD,"/mZ=",mZ_MuSTARMD,"/TZ=",TZ_MuSTARMD,"/",pname_MuSTARMD,"/nEX=",numEX_MuSTARMD,"/fq=",fq_MuSTARMD,"/pmf/pmf_TAA=",TAA_MuSTARMD,"_TCG_",TCG_MuSTARMD,"_TZ_",TZ_MuSTARMD,"_",width_MuSTARMD,"_",KZAAo,"_",KZCGo,"_",AACG,"_4_2012-08-21",sep="") 
felwFillConMapwrangewy(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,mar=c(5,5,2,2))

felwFillConBar(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

eof

Rscript ${paraRfil}
