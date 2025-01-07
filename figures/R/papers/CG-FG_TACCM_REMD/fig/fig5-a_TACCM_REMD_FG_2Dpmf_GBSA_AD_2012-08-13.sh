
dir=~/papers/CG-FG_TACCM_REMD/

parafile=~/calspa/InV2InW/AD/para/para_e_CG-FG_NH_2012-07-27_wrefd0.1_TZ=750_fq=10ps_2_GBSA_2012-08-06_ff99SB.sh

source ${parafile}

paraRfil=${dir}/fig/para_fig5-a_TACCM_REMD_FG_2Dpmf_GBSA_AD_2012-08-13.R

cat <<EOF > ${paraRfil}

TAA<-"${TAA}"
TCG<-"${TCG}"

TZ<-"${TZ}"

tau<- "1.0"

mZ <- "50000.00"

TZ <- "${TZs[1]}"

pname <- "${pname}"

KZAAo <- "${KZAAo[2]}"

KZCGo <- "${KZCGo[2]}"

nEX <- "${numEX}"

AACG <- "${MODE[2]}"

fq <- "${TLbase}"

level <- seq(0,20,2)

width <- "0.3"

name.title <- " "

WVflag <- "InW"

name.out <- "~/papers/CG-FG_TACCM_REMD/tiff/fig5-a_TACCM_REMD_FG_2Dpmf_GBSA_AD_2012-08-13"

source("~/papers/CG-FG_TACCM_REMD/graph/graph_TACCM_REMD_2Dpmf_GBSA_AD_2012-08-13.R")

EOF

Rscript ${paraRfil}; echo ${paraRfil}
