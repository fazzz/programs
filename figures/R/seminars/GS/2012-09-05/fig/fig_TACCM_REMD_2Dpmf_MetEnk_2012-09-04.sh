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

dir=~/calspa/TACCM_CGAAREMD/MetEnk

parafile=${dir}/para/para_e_CG-FG_NH_2CG21_KAA=1000_TZ=3000_TCG=370_freq=10ps_2012-08-30.sh
origAA=0
origCG1=0
origCG2=0
heightAA=16
heightCG1=16
heightCG2=16
dwidth=1
width1=0.005

source ${parafile}
TL=`expr  ${TLbase}`
TS=`expr ${TLbase} \* ${numEX}`

tiffheight=$( echo ${tiffwidth} \* 320 | bc )
tiffheight=$( echo ${tiffheight} / 390 | bc )

echo width=${tiffwidth}
echo height=${tiffheight}

AACGflag=( dummy AA CG1 CG2 )
orig=( dummy )
orig[1]=${origAA} 
orig[2]=${origCG1}
orig[3]=${origCG2}

height=( dummy )
height[1]=${heightAA} 
height[2]=${heightCG1}
height[3]=${heightCG2}

dwidth[1]=${dwidth}
dwidth[2]=${dwidth}
dwidth[3]=${dwidth}

paraRfil=( dummy ${dir}/R/para_e_2CG1FG_CMD_NH_${pname}_TZ=${TZ}_freq_${freq}ps_AA_wTS_2012-08-30.R ${dir}/R/para_e_2CG1FG_CMD_NH_${pname}_TZ=${TZ}_freq_${freq}ps_CG1_wTS_2012-08-30.R ${dir}/R/para_e_2CG1FG_CMD_NH_${pname}_TZ=${TZ}_freq_${freq}ps_CG2_wTS_2012-08-30.R )

for i in `seq 1 3`; do
    cat <<EOF > ${paraRfil[$i]}
TAA<-"${TAA}"
TCG1<-"${TCG1}"
TCG2<-"${TCG2}"
TCG<-"${TCG}"

TZs <- c( "${TZs[1]}" )
nTZs <- length(TZs)

tau<-c(  "1.0" )
ntau<-length(tau)

ep<-"${ep}"
nb<-"${nb}"
cutoff<-"${cutoff}"

mZ<-c(  "50000.00" )
nmZ<-length(mZ)

nKZAA<-${nKZAA}
nKZCG1<-${nKZCG1}
nKZCG2<-${nKZCG2}
numRE<-${numRE}

pname<-"${pname}"

numEX<-"${numEX}"

AACG <- "${AACGflag[$i]}"
width <- "${width1}"

KZAAo <- "${KZAAo[$i]}"
KZCG1o <- "${KZCG1o[$i]}"
KZCG2o <- "${KZCG2o[$i]}"

level <- seq(${orig[$i]},${height[$i]},${dwidth[$i]})

state <- "1FG2CG"

TS <- "${TS}"

TLbase <- "${TLbase}"

name.title <- paste(AACG,sep="")

tiffwidth <- ${tiffwidth}

tiffheight <- ${tiffheight}

name.out <- "~/seminars/GS/2012-09-05/tiff/fig_pmf2D_TACCM_REMD_MetEnk_e_CG-FG_NH_2012-09-04_${AACGflag[$i]}"

source("~/seminars/GS/2012-09-05/graph/graph_TACCM_REMD_2Dpmf_MetEnk_2012-09-04.R")
EOF
    Rscript ${paraRfil[$i]}; echo ${paraRfil[$i]}
done

