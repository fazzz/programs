#!/bin/zsh

dir=~/Report/2013-08/

parafile1=~/calspa/refcalc/REMD/AD/para/2013-08-17/para_s_REMD_vac_ff99SB_2013-08-17_TZ=300-600.sh
parafile2=~/calspa/refcalc/REMD/AD/para/2013-08-17/para_s_REMD_vac_ff99SB_2013-08-17_TZ=300-800.sh
parafile3=~/calspa/refcalc/REMD/AD/para/2013-08-17/para_s_REMD_vac_ff99SB_2013-08-17_TZ=300-1000.sh
parafile4=~/calspa/refcalc/REMD/AD/para/2013-08-17/para_s_REMD_vac_ff99SB_2013-08-17_TZ=300-1200.sh

paraRfil=${dir}/graph/fig_a-e_40ns-REMD_2Dpmf_AD_2013-08-17.R

cat << eof > ${paraRfil}
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

file.name <- "~/Report/2013-08/eps/fig_a-e_40ns-REMD_2Dpmf_AD_2013-08-17.eps"

postscript(file.name,width=4.0,height=9.0,horizontal=FALSE,onefile=FALSE,paper="special")

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

source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange_worwoaxis.R")
source("~/papers/CG-FG_TACCM_REMD/FillConBar.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap_wrange_worwoaxis_wpoints.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConBar_wmai.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMapbox_wpoints_2.R")

nf <- layout(matrix(c(1,5,2,5,3,5,4,5),4,2,byrow=TRUE),c(2.9,1.1),c(2.5,2.0,2.0,2.8))

par(cex.axis=1.0)
par(cex.lab=1.0)

label.x <- ""
label.y <- ""

title=name.title
eof

source ${parafile1}
cat << eof >> ${paraRfil}
pname_TREMD <-"${pname}"
numEX_TREMD <-"${numEX}"
TLbase_TREMD <-"${TLbase}"
ff_TREMD <-"${ff}"

dir0 <- "~/calspa/refcalc/REMD/AD"
dirbase <- paste(dir0,"/s_REVAC_2012-07-23_",ff_TREMD,sep="")
name <- paste(dirbase,"/",pname_TREMD,"/nEX=",numEX_TREMD,"/freq=",TLbase_TREMD,"ps","/pmf/pmf_py3",sep="")

cat(name)
	
felwFillConMapwrangeworwoaxis(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,xaflag="f",mai=c(0.0,0.8,0.5,0.1),norm="T")
eof

source ${parafile2}
cat << eof >> ${paraRfil}
pname_TREMD <-"${pname}"

numEX_TREMD <-"${numEX}"

TLbase_TREMD <-"${TLbase}"

ff_TREMD <-"${ff}"

dir0 <- "~/calspa/refcalc/REMD/AD"
dirbase <- paste(dir0,"/s_REVAC_2012-07-23_",ff_TREMD,sep="")
name <- paste(dirbase,"/",pname_TREMD,"/nEX=",numEX_TREMD,"/freq=",TLbase_TREMD,"ps","/pmf/pmf_py3",sep="")

cat(name)
	
felwFillConMapwrangeworwoaxis(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,xaflag="f",mai=c(0.0,0.8,0.0,0.1),norm="T")
eof

source ${parafile3}
cat << eof >> ${paraRfil}
pname_TREMD <-"${pname}"

numEX_TREMD <-"${numEX}"

TLbase_TREMD <-"${TLbase}"

ff_TREMD <-"${ff}"

dir0 <- "~/calspa/refcalc/REMD/AD"
dirbase <- paste(dir0,"/s_REVAC_2012-07-23_",ff_TREMD,sep="")
name <- paste(dirbase,"/",pname_TREMD,"/nEX=",numEX_TREMD,"/freq=",TLbase_TREMD,"ps","/pmf/pmf_py3",sep="")

cat(name)
	
felwFillConMapwrangeworwoaxis(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,xaflag="f",mai=c(0.0,0.8,0.0,0.1),norm="T")
eof

source ${parafile4}
cat << eof >> ${paraRfil}
pname_TREMD <-"${pname}"

numEX_TREMD <-"${numEX}"

TLbase_TREMD <-"${TLbase}"

ff_TREMD <-"${ff}"

dir0 <- "~/calspa/refcalc/REMD/AD"
dirbase <- paste(dir0,"/s_REVAC_2012-07-23_",ff_TREMD,sep="")
name <- paste(dirbase,"/",pname_TREMD,"/nEX=",numEX_TREMD,"/freq=",TLbase_TREMD,"ps","/pmf/pmf_py3",sep="")

cat(name)
	
felwFillConMapwrangeworwoaxis(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,mai=c(0.8,0.8,0.0,0.1),norm="T")
eof

cat <<eof >> ${paraRfil}

felwFillConBarwmai(name,label.x=label.x,label.y=label.y,level=level,mai=c(7.2,0.5,0.5,0.2))

eof

Rscript ${paraRfil}
