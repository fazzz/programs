#!~/bin/sh

parafile=~/calspa/TACCM/AD/para/para_e_ff99SB_NH_2012-07-31_TZ=750.sh

dirbase=~/papers/CG-FG_TACCM_REMD

source ${parafile}

paraRfil=${dirbase}/fig/para_pmf2D_s_TAMD_VAC_2012-08-08.R

i=1

cat <<EOF > ${paraRfil}
T<-"${T}"
TB<-"${TB[$i]}"

tau<-"${tau}"

KZ<-"${KZ}"

mZ<-"${mZ}"

width <- "0.3"

level <- seq(0,20,2)

name.out <- "${dirbase}/tiff/fig4-e_TAMD_2Dpmf_AD_2012-08-08"

ff <- "${ff}"

source("~/papers/CG-FG_TACCM_REMD/graph/graph_TAMD_2Dpmf_AD_2012-08-08.R")
EOF
Rscript ${paraRfil}; echo ${paraRfil}
