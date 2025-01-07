#!/bin/zsh

dir=~/papers/CG-FG_TACCM_REMD/

parafile=~/calspa/TACCM_CGAAREMD/AD/para/para_e_CG-FG_NH_2012-07-27_wrefd0.1_TZ=750_fq=10ps_99SB.sh 

source ${parafile}

paraRfil=${dir}/fig/para_fig3_a-b-c_MuSTARMD_AD_2Dpmf_2012-10-03.R

tiffwidth=650
tiffheight=200

cat << eof > ${paraRfil}
TAA<-"${TAA}"
TCG<-"${TCG}"

TZ <-"${TZs[1]}"

tau<-"${tau[1]}"

mZ <-"${mZ[1]}"  

pname<-"${pname}"

numEX<-"${numEX}"

fq<-"${TLbase}"

TLbase<-"${TLbase}"

width <- "0.3"

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

file.name <- "~/papers/CG-FG_TACCM_REMD/tiff/fig3_a-b-c_MuSTARMD_AD_2Dpmf_2012-10-03.tiff"
tiff(file.name,width=${tiffwidth},height=${tiffheight})

par(oma = c(0,0,0,0) )
source("~/papers/CG-FG_TACCM_REMD/FillConMap.R")
source("~/papers/CG-FG_TACCM_REMD/FillConBar.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConBar.R")
						
nf <- layout(matrix(c(1,2,3,4),1,4,byrow=TRUE),c(3,3,3,1),c(1))

par(mar =c(5,5,2,2) )
par(cex.axis=2.0)
par(cex.lab=2.0)
label.x <- expression(paste(phi))
label.y <- expression(paste(psi))
dir0 <- "/home/yamamori/calspa/TACCM_CGAAREMD/AD"
dirbase <- paste(dir0,"/e_CG-FG_NH_2012-07-27",sep="")

title=name.title
	
AACG <- "AA"
KZAAo <- "1000"
KZCGo <- "0"   
name <- paste(dirbase,"/tau=",tau,"/mZ=",mZ,"/TZ=",TZ,"/",pname,"/nEX=",numEX,"/fq=",fq,"/pmf/pmf_TAA=",TAA,"_TCG_",TCG,"_TZ_",TZ,"_",width,"_",KZAAo,"_",KZCGo,"_",AACG,"_4_2012-08-21",sep="") 
felwFillConMap(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

AACG <- "CG"
KZAAo <- "0"
KZCGo <- "1000"
name <- paste(dirbase,"/tau=",tau,"/mZ=",mZ,"/TZ=",TZ,"/",pname,"/nEX=",numEX,"/fq=",fq,"/pmf/pmf_TAA=",TAA,"_TCG_",TCG,"_TZ_",TZ,"_",width,"_",KZAAo,"_",KZCGo,"_",AACG,"_4_2012-08-21",sep="")
felwFillConMap(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")
														   
AACG <- "AA"
KZAAo <- "1000"
KZCGo <- "0"   
name <- paste(dirbase,"/tau=",tau,"/mZ=",mZ,"/TZ=",TZ,"/",pname,"/nEX=",numEX,"/fq=",fq,"/pmf/pmf_TAA=",TAA,"_TCG_",TCG,"_TZ_",TZ,"_",width,"_",KZAAo,"_",KZCGo,"_",AACG,"_4_2012-08-21_0",sep="")
felwFillConMap(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

felwFillConBar(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

eof

Rscript ${paraRfil}

