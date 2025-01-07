#!/bin/zsh

dir=~/Report/2013-08/

paraRfil=${dir}/graph/fig_80ns-REUS_2Dpmf_AD_2013-08-17.R

tiffwidth=720
tiffheight=300

parafile_REUS=~/calspa/refcalc/REUS/AD/para/para_s_REUSMD_vac_ff99SB_Nbin=4x4_K=1_10ns_freq=10ps_2013-08-12.sh
source ${parafile_REUS}
cat << eof > ${paraRfil}
fact.x=1
fact.y=1

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

file.name <- "~/Report/2013-08/eps/fig_80ns-REUS_2Dpmf_AD_2013-08-17.eps"

postscript(file.name,width=3.6,height=3.0,horizontal=FALSE,onefile=FALSE,paper="special")

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

nf <- layout(matrix(c(1,2),1,2,byrow=TRUE),c(3,1),c(1))

par(cex.axis=1.0)
par(cex.lab=1.0)

label.x <- expression(paste(phi,"(radian/",pi,")"))
label.y <- expression(paste(,psi,"(radian/",pi,")"))

title=name.title

pname_REUS <-"${pname}"

numEX_REUS <-"${numEX}"

TLbase_REUS <-"${TLbase}"

ff_REUS <-"${ff}"

dir0 <- "~/calspa/refcalc/REUS/AD"
dirbase <- paste(dir0,"/s_REUSVAC_2013-08-12_",ff_REUS,sep="")

name <- paste(dirbase,"/",pname_REUS,"/nEX_",numEX_REUS,"/freq_",TLbase_REUS,"ps","/pmf/pmf_UmbMD_vac_nbx=40_nby=40.txt_2",sep="")

cat(name)
	
felwFillConMapwrangeworwoaxis(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,mai=c(0.8,0.8,0.5,0.1),norm="T")

felwFillConBarwmai(name,label.x=label.x,label.y=label.y,level=level,mai=c(0.8,0.5,0.5,0.2))

eof

Rscript ${paraRfil}
