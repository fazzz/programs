#!/bin/zsh

dir=~/papers/CG-FG_TACCM_REMD/

paraRfil=${dir}/fig/para_fig4_a-b-c-d_comp_MuSTARMD-other_methods_AD_2Dpmf_2012-10-03.R

tiffwidth=700
tiffheight=550

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

parafile_TREMD=~/calspa/refcalc/REMD/AD/para/para_s_REMD_vac_ff99SB_2012-07-24.sh
source ${parafile_TREMD}
cat <<eof >> ${paraRfil}
pname_TREMD <-"${pname}"

numEX_TREMD <-"${numEX}"

TLbase_TREMD <-"${TLbase}"

ff_TREMD <-"${ff}"

eof

parafile_TAMD=~/calspa/TACCM/AD/para/para_e_ff99SB_NH_2012-07-31_TZ=750.sh
source ${parafile_TAMD}
i=1
cat <<eof >> ${paraRfil}
T_TAMD<-"${T}"

TB_TAMD<-"${TB[$i]}"

tau_TAMD<-"${tau}"

KZ_TAMD<-"${KZ}"

mZ_TAMD<-"${mZ}"

width_TAMD <- "0.3"

ff_TAMD <- "${ff}"

eof

parafile_CMD=~/calspa/refcalc/CMD/AD/para/para_s_CMD_T300_ff99SB_2012-07-24.sh
source ${parafile_CMD}
cat <<eof >> ${paraRfil}
pname_CMD<-"${pname}"

T_CMD<-"${T}"

ff_CMD<-"${ff}"

width_CMD <- "0.2"

eof

parafileUmbSam=~/calspa/refcalc/UmbSam/AD/para/para_s_vac_ff99SB_K=10_2012-11-12.sh
source ${parafileUmbSam}
TLns=`expr ${TLbase} / 1000`
cat <<eof >> ${paraRfil}
pname_UMB<-"${pname}"

ff_UMB<-"${ff}"

TLns<-"${TLns}"

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

file.name <- "~/papers/CG-FG_TACCM_REMD/eps/fig3_a-b-c-d_comp_MuSTARMD-other_methods_AD_2Dpmf_2013-03-12.tiff"
tiff(file.name,width=${tiffwidth},height=${tiffheight})
#postscript(file.name,width=5.6,height=4.5,horizontal=FALSE,onefile=FALSE,paper="special")

xrange <- c(-3,3,6)
yrange <- c(-3,3,6)

par(oma = c(0,0,0,0) )
source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange.R")
source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange_wy.R")
source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange_wox_wy.R")
source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange_woxyaxis.R")

source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap_wrange.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap_wrange_wy.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap_wrange_woxwy.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap_wrange_woxyaxis.R")

source("~/papers/CG-FG_TACCM_REMD/FillConBar.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConBar.R")

nf <- layout(matrix(c(1,2,3,4,5,6),2,3,byrow=TRUE),c(3.8,3,3.8),c(1,1.1))

#par(mar =c(5,5,2,2) )

par(cex.axis=1.0)
par(cex.lab=2.0)
label.x <- expression(paste(phi))
label.y <- expression(paste(psi))

title=name.title
	
dir0 <- "/home/yamamori/calspa/TACCM_CGAAREMD/AD"
dirbase <- paste(dir0,"/e_CG-FG_NH_2012-07-27",sep="")
AACG <- "CG"
KZAAo <- "0"
KZCGo <- "1000"
name <- paste(dirbase,"/tau=",tau_MuSTARMD,"/mZ=",mZ_MuSTARMD,"/TZ=",TZ_MuSTARMD,"/",pname_MuSTARMD,"/nEX=",numEX_MuSTARMD,"/fq=",fq_MuSTARMD,"/pmf/pmf_TAA=",TAA_MuSTARMD,"_TCG_",TCG_MuSTARMD,"_TZ_",TZ_MuSTARMD,"_",width_MuSTARMD,"_",KZAAo,"_",KZCGo,"_",AACG,"_4_2012-08-21",sep="") 

felwFillConMapwrangewy(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,mar=c(5,5,1,0))

dir0 <- "/home/yamamori/calspa/TACCM/AD"
dirbase <- paste(dir0,"/e_TACCM_NH_2012-07-31_",ff_TAMD,sep="")
name <- paste(dirbase,"/tau=",tau_TAMD,"/TB=",TB_TAMD,"/KZ=",KZ_TAMD,"/mZ=",mZ_TAMD,"/pmf/pmf_T=",T_TAMD,".Zhist",sep="")

felwFillConMapwrange(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,mar=c(5,0,1,0),plot.axis="yes")

dir0 <- "~/calspa/refcalc/REMD/AD"
dirbase <- paste(dir0,"/s_REVAC_2012-07-23_",ff_TREMD,sep="")
name <- paste(dirbase,"/",pname_TREMD,"/nEX=",numEX_TREMD,"/freq=",TLbase_TREMD,"ps","/pmf/pmf_ADv",sep="")

felwFillConMapwrange(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,mar=c(5,0,1,5),plot.axis="yes")

dir0 <- "~/calspa/refcalc/CMD/AD"
dirbase <- paste(dir0,"/s_CVAC_2012-08-08_",ff_CMD,sep="")
name <- paste(dirbase,"/anl/pmf_ADv_T",T_CMD,"_",width_CMD,sep="")

felwFillConMapwrangewy(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,mar=c(5,5,2,0),plot.axis="yes")

dir0 <- "~/calspa/refcalc/UmbSam/AD"
dirbase <- paste(dir0,"/s_UmbSam_vac_2012-11-12_",ff_UMB,sep="")
name <- paste(dirbase,"/",pname_UMB,"/pmf/pmf_UmbMD_vac_",TLns,"ns.txt",sep="")

#dir0 <- "~/calspa/refcalc/CMD/AD"
#dirbase <- paste(dir0,"/s_CVAC_2012-08-08_",ff_CMD,sep="")
#name <- paste(dirbase,"/anl/pmf_ADv_T",T_CMD,"_",width_CMD,sep="")

felwFillConMapwrange(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,mar=c(5,0,2,0),plot.axis="yes")

par(mar=c(5,20,2,8))

felwFillConBar(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

eof

Rscript ${paraRfil}

