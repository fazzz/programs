#!/bin/zsh

dir=~/papers/CG-FG_TACCM_REMD/

paraRfil=${dir}/graph/fig7_a-b-c_comp_MuSTARMD-TREMD_CMD_MetEnk_2Dpmf_2013-02-28.R

parafile=~/calspa/TACCM_CGAAREMD/MetEnk/para/2013-05-31/para_e_CG-FG_NH_2CG_KAA=3000_TZ=3000_freq=1ps_2013-06-01.sh

source ${parafile}

tiffwidth=800
tiffheight=230

cat << eof > ${paraRfil}
fact.x <- 1
fact.y <- 1

#max <- "0.4"
max <- "0.3"

TAA<-"${TAA}"
TCG1<-"${TCG1}"
TCG2<-"${TCG2}"

TZs <-"${TZs[1]}"
tau<-"${tau[1]}"

ep<-"${ep}"
nb<-"${nb}"
cutoff<-"${cutoff}"

mZ<-"${mZ[1]}"

nKZAA<-${nKZAA}
nKZCG1<-${nKZCG1}
nKZCG2<-${nKZCG2}
numRE<-${numRE}

pname<-"${pname}"

freq<-"${TLbase}"

AACG <- "AA"
width <- "0.01"

KZAAo <- "${KZAAo[1]}"
KZCG1o <- "${KZCG1o[1]}"
KZCG2o <- "${KZCG2o[1]}"

level <- seq(0,10,1)

state <- "1FG2CG"

TS <- "${numEX}"

TLbase <- "${TLbase}"

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

file.name <- "~/papers/CG-FG_TACCM_REMD/eps/fig7_a-b-c_comp_MuSTARMD-TREMD_CMD_MetEnk_2Dpmf_2013-05-31.eps"
#tiff(file.name,width=${tiffwidth},height=${tiffheight})
postscript(file.name,width=3.2,height=5.8,horizontal=FALSE,onefile=FALSE,paper="special")

par(oma = c(0,0,0,0) )
source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange_worwoaxis.R")
source("~/papers/CG-FG_TACCM_REMD/FillConBar.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap_wrange_worwoaxis.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConBar_wmai.R")

#xrange <- c(0,0.3,3)

#nf <- layout(matrix(c(1,2,3,4),1,4,byrow=TRUE),c(3.5,3,3,1.0),c(1))
nf <- layout(matrix(c(1,4,2,4,3,4),3,2,byrow=TRUE),c(3,1),c(20,15,23))

par(cex.axis=1.2)
par(cex.lab=1.2)

label.x <- "d1"
label.y <- "d2"

name.title <- NULL

title=name.title
	
dir0 <- "/home/yamamori/calspa/TACCM_CGAAREMD/MetEnk/"
dirbase <- paste(dir0,"/e_CG-FG_CMD_NH_2012-06-04/",state,sep="")

nTZs<-1

#name <- paste(dirbase,"/tau=",tau,"/mZ=",mZ,"/ep=",ep,"/cutoff=",cutoff,"/TZ=",TZs,"/",pname,"/freq=",freq,"/pmf/pmf_TAA=",TAA,"_TCG1_",TCG1,"_TCG2_",TCG2,"_TZ_",TZs,"_",width,"_",KZAAo,"_",KZCG1o,"_",KZCG2o,"_",AACG,"_bo10000ps",sep="")

name <- paste(dirbase,"/tau=",tau,"/mZ=",mZ,"/ep=",ep,"/cutoff=",cutoff,"/TZ=",TZs,"/",pname,"/freq=",freq,"/pmf/pmf_pymbar_AA_41_2",sep="")

cat(name)

#yrange <- c(0.1,0.3,2)
#yrange.axis <- c(0.1,0.3,2)

label.y <- " "

#xrange <- c(0,0.4,4)
#xrange.axis <- c(0,0.4,4)
#yrange <- c(0,0.4,4)
#yrange.axis <- c(0.1,0.4,3)

xrange <- c(0,0.3,3)       
xrange.axis <- c(0,0.3,3)  
yrange <- c(0,0.3,3)       
yrange.axis <- c(0.1,0.3,2)

felwFillConMapwrangeworwoaxis(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,mai=c(0.0,0.8,0.5,0.1),xaflag="f")

text(0.03,0.27,"(a)",cex=1.6)

#xrange <- c(0.05,0.3,5)
#xrange.axis <- c(0.05,0.3,5)
#yrange <- c(0.05,0.3,5)
#yrange.axis <- c(0.05,0.3,5)

#name <- "/home/yamamori/calspa/TACCM_CGAAREMD/MetEnk//e_CG-FG_NH_2012-05-14/FG_Amber/tau=1.0/anl/pmf_MetEnk_AA_T=300_fLM1_0.005.txt"

#name <- "/home/yamamori/calspa/refcalc/CMD/MetEnk/s_CVAC_2013-03-13_LM1/anl/pmf_MetEnk_T=300_fLM1_0.005.txt"

name <- "/home/yamamori/calspa/refcalc/CMD/MetEnk/s_CVAC_2013-03-13_LM1/anl/pmf_MetEnk_T=300_fLM1_0.01.txt"

#name <- "/home/yamamori/calspa/refcalc/CMD/MetEnk/s_CVAC_2013-03-13_LM1/anl/pmf_MetEnk_T=300_fLM1_0.01_0.0-0.4.txt"

cat(name)

#yrange <- c(0.1,0.3,2)
#yrange.axis <- c(0.1,0.3,2)

label.y <- "d2"

#xrange <- c(0,0.4,4)
#xrange.axis <- c(0,0.4,4)
#yrange <- c(0,0.4,4)
#yrange.axis <- c(0.1,0.4,3)

xrange <- c(0,0.3,3)
xrange.axis <- c(0,0.3,3)
yrange <- c(0,0.3,3)
yrange.axis <- c(0.1,0.3,2)

felwFillConMapwrangeworwoaxis(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,mai=c(0.0,0.8,0.0,0.1),xaflag="f")

text(0.03,0.27,"(b)",cex=1.6)

#xrange <- c(0.05,0.3,5)
#xrange.axis <- c(0.05,0.3,5)
#yrange <- c(0.05,0.3,5)
#yrange.axis <- c(0.05,0.3,5)

#name <- "/home/yamamori/calspa/TACCM_CGAAREMD/MetEnk//e_CG-FG_NH_2012-05-14/FG_Amber/tau=1.0/anl/pmf_MetEnk_AA_T=300_fLM2_0.005.txt"

#name <- "/home/yamamori/calspa/refcalc/CMD/MetEnk/s_CVAC_2013-03-13_LM2/anl/pmf_MetEnk_T=300_fLM2_0.005.txt"

name <- "/home/yamamori/calspa/refcalc/CMD/MetEnk/s_CVAC_2013-03-13_LM2/anl/pmf_MetEnk_T=300_fLM2_0.01.txt"

#name <- "/home/yamamori/calspa/refcalc/CMD/MetEnk/s_CVAC_2013-03-13_LM2/anl/pmf_MetEnk_T=300_fLM2_0.01_0.0-0.4.txt"

cat(name)

label.y <- " "

#xrange <- c(0,0.4,4)
#xrange.axis <- c(0,0.4,4)
#yrange <- c(0,0.4,4)
#yrange.axis <- c(0,0.4,4)

xrange <- c(0.0,0.3,3)
xrange.axis <- c(0.0,0.3,3)
yrange <- c(0.0,0.3,3)
yrange.axis <- c(0.0,0.3,3)

felwFillConMapwrangeworwoaxis(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,mai=c(0.8,0.8,0.0,0.1))

par(cex.axis=1.0)

text(0.03,0.27,"(c)",cex=1.6)

#felwFillConBarwmai(name,label.x=label.x,label.y=label.y,level=level,mai=c(0.8,0.5,0.5,0.2))
felwFillConBarwmai(name,label.x=label.x,label.y=label.y,level=level,mai=c(0.8,0.5,0.5,0.5))

eof

Rscript ${paraRfil}
