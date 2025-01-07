#!~/bin/sh

parafile=~/calspa/refcalc/REMD/AD/para/para_s_REMD_vac_ff99SB_2012-07-24.sh

dirbase=~/papers/CG-FG_TACCM_REMD

source ${parafile}

paraRfil=~/papers/CG-FG_TACCM_REMD/fig/para_fig4-d_REMD_2Dpmf_AD_2012-08-08.R

cat <<EOF > ${paraRfil}
pname<-"${pname}"

numEX<-"${numEX}"

TLbase<-"${TLbase}"

level <- seq(0,20,2)

ff<-"${ff}"

name.out <- paste("${dirbase}/tiff/fig4-d_REMD_2Dpmf_AD_2012-08-08",sep="")

source("${dirbase}/graph/graph_REMD_2Dpmf_AD_2012-08-08.R")
EOF
Rscript ${paraRfil}; echo ${paraRfil}


