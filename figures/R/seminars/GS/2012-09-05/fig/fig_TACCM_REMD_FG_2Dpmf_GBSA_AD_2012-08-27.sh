#!~/bin/sh

opt=(dummy tiffwidth )
nopt=${#opt[*]}
if [ $# -le `expr ${nopt} - 2`  ]; then
    echo "USAGE " $0  ${opt[*]:1:${nopt}}
    echo $*
    exit
fi

num=1
while [ $num -le `expr ${nopt} - 1` ]; do
    eval ${opt[$num]}=$1
    shift 1
    num=`expr $num + 1`
done

dir=~/seminars/GS/2012-09-05
#parafile=~/calspa/InV2InW/AD/para/para_e_CG-FG_NH_2012-07-27_wrefd0.1_TZ=750_fq=10ps_2_GBSA_2012-08-06_ff99SB.sh
parafile=~/calspa/InV2InW/AD/para/para_e_CG-FG_NH_2012-07-27_wrefd0.1_TZ=750_fq=10ps_2_GBSA_2012-08-06.sh

source ${parafile}

tiffheight=$( echo ${tiffwidth} \* 320 | bc )
tiffheight=$( echo ${tiffheight} / 390 | bc )

echo width=${tiffwidth}
echo height=${tiffheight}

paraRfil=${dir}/fig/para_fig5-a_TACCM_REMD_FG_2Dpmf_GBSA_AD_2012-08-27.R

cat <<EOF > ${paraRfil}

TAA<-"${TAA}"
TCG<-"${TCG}"

TZ<-"${TZ}"

tau<- "1.0"

mZ <- "50000.00"

TZ <- "${TZs[1]}"

pname <- "${pname}"

KZAAo <- "${KZAAo[1]}"

KZCGo <- "${KZCGo[1]}"

nEX <- "${numEX}"

AACG <- "AA"

fq <- "${TLbase}"

level <- seq(0,20,2)

width <- "0.3"

name.title <- " "

tiffwidth <- ${tiffwidth}

tiffheight <- ${tiffheight}

WVflag <- "InW"

name.out <- "${dir}/tiff/fig_TACCM_REMD_FG_2Dpmf_GBSA_AD_2012-08-27"

source("${dir}/graph/graph_TACCM_REMD_2Dpmf_GBSA_AD_2012-08-27.R")

EOF

Rscript ${paraRfil}; echo ${paraRfil}
