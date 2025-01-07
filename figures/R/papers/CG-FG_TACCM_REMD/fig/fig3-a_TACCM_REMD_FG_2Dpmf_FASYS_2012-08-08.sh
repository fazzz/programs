#!~/bin/sh

parafile=~/calspa/TACCM_CGAAREMD/FiveAtomSysf/para/para_e_CG-FG_NH_2012-08-14_T1_TZ=1000_FG7-CG4_40x8ns.sh

dir=~/calspa/TACCM_CGAAREMD/FiveAtomSysf

source ${parafile}

AACGflag=CG
height=10

dwidth=1

paraRfil=${dir}/R/para_fig3-a_TACCM_REMD_FG_2Dpmf_FASYS_2012-08-08.R

cat <<EOF > ${paraRfil}
TAA<-"${TAA}"
TCG<-"${TCG}"

TZs <- "${TZs[1]}"
nTZs <- length(TZs)

tau<-c(  "1.0" )
ntau<-length(tau)

mZ<-c(  "50000.00" )
nmZ<-length(mZ)

nKZAA<-${nKZAA}
nKZCG<-${nKZCG}
numRE<-${numRE}

pname<-"${pname}"

numEX<-"${numEX}"

fq<-"${TLbase}"

AACG <- "AA"
width <- "0.3"

KZAAo <- "1000"
KZCGo <- "0"

level <- seq(0,${height},${dwidth})

name.title <- NULL

name.out <- "~/papers/CG-FG_TACCM_REMD/tiff/fig3-a_TACCM_REMD_FG_2Dpmf_FASYS_2012-08-08"

source("~/papers/CG-FG_TACCM_REMD/graph/graph_TACCM_REMD_2Dpmf_FASYS_2012-08-08.R")
EOF

Rscript ${paraRfil}; echo ${paraRfil}
