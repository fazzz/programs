#!~/bin/sh

parafile=~/calspa/refcalc/REMD/AD/para/para_s_REMD_GBSA_ff99SB_2012-07-24.sh

dirbase=~/papers/CG-FG_TACCM_REMD

source ${parafile}

paraRfil=~/papers/CG-FG_TACCM_REMD/fig/para_fig5-b_REMD_2Dpmf_GBSA_AD_2012-08-13.R

cat <<EOF > ${paraRfil}
pname<-"${pname}"

numEX<-"${numEX}"

TLbase<-"${TLbase}"

level <- seq(0,20,2)

ff<-"${ff}"

name.out <- paste("${dirbase}/tiff/fig5-b_REMD_2Dpmf_GBSA_AD_2012-08-13",sep="")

source("${dirbase}/graph/graph_REMD_2Dpmf_GBSA_AD_2012-08-13.R")
EOF

Rscript ${paraRfil}; echo ${paraRfil}


