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

parafile=~/calspa/TACCM_CGAAREMD/FiveAtomSysf/para/para_e_CG-FG_NH_2012-08-14_T1_TZ=1000_FG7-CG4_40x8ns.sh

dir=~/seminars/GS/2012-09-05/

source ${parafile}

tiffheight=$( echo ${tiffwidth} \* 320 | bc )
tiffheight=$( echo ${tiffheight} / 390 | bc )

echo width=${tiffwidth}
echo height=${tiffheight}

AACGflag=CG
height=10

dwidth=1

paraRfil=${dir}/fig/para_fig_TACCM_REMD_FG_2Dpmf_FASYS_2012-08-27.R

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

name.out <- "~/seminars/GS/2012-09-05/tiff/fig_TACCM_REMD_FG_2Dpmf_FASYS_2012-08-27"

tiffwidth <- ${tiffwidth}

tiffheight <- ${tiffheight}

source("~/seminars/GS/2012-09-05/graph/graph_TACCM_REMD_2Dpmf_FASYS_2012-08-27.R")
EOF

Rscript ${paraRfil}; echo ${paraRfil}
