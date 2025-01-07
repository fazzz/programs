#!~/bin/sh

opt=(dummy size )
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

dir=~/calspa/TACCM_CGAAREMD/AD

heightAA=10
heightCG=10
dwidth=1
width1=0.3 
numuene=4
parafile=~/calspa/TACCM_CGAAREMD/AD/para/2012-09-17/para_e_CG-FG_NR=4_wrefd0.1_TZ=350_fq=10ps_99SB_KZmax=1000.sh

source ${parafile}

TZ=${TZs[1]}

AACGflag=( dummy AA CG )
height=( dummy )
height[1]=${heightAA} 
height[2]=${heightCG}

dwidth[1]=${dwidth}
dwidth[2]=${dwidth}

paraRfil=~/Report/2012-09/graph/para_2Dpmf_TACCM_CGAA_${pname}_TZ=${TZ}_2012-09-25.R

cat <<EOF > ${paraRfil}
TAA<-"${TAA}"
TCG<-"${TCG}"

TZs <- c( "${TZ}" )
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

TLbase<-"${TLbase}"

AACG <- "${AACGflag[2]}"
width <- "${width1}"

KZAAo <- "${KZAAo[2]}"
KZCGo <- "${KZCGo[2]}"

numuene <- "${numuene}"

level <- seq(0,${height[2]},${dwidth[2]})

name.title <- paste(AACG,sep="")

name.out <- "~/Report/2012-09/tiff/pmf_2D_TACCM_CGAA_pymbar_${pname}_TZ=${TZ}_2012-09-25"

size <- ${size}

source("~/Report/2012-09/graph/graph_2Dpmf_TACCM_CGAA_pymbar_2012-09-25.R")
EOF
Rscript ${paraRfil}; echo ${paraRfil}
