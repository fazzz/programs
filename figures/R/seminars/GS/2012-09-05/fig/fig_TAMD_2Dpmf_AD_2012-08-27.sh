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

parafile=~/calspa/TACCM/AD/para/para_e_ff99SB_NH_2012-07-31_TZ=750.sh

dirbase=~/seminars/GS/2012-09-05

source ${parafile}

tiffheight=$( echo ${tiffwidth} \* 320 | bc )
tiffheight=$( echo ${tiffheight} / 390 | bc )

echo width=${tiffwidth}
echo height=${tiffheight}

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

name.out <- "${dirbase}/tiff/fig_TAMD_2Dpmf_AD_2012-08-27"

ff <- "${ff}"

tiffwidth <- ${tiffwidth}

tiffheight <- ${tiffheight}

source("${dirbase}/graph/graph_TAMD_2Dpmf_AD_2012-08-27.R")
EOF
Rscript ${paraRfil}; echo ${paraRfil}
