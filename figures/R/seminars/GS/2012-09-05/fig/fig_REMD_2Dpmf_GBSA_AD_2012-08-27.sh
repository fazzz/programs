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
#parafile=~/calspa/refcalc/REMD/AD/para/para_s_REMD_GBSA_ff99SB_2012-07-24.sh
parafile=~/calspa/refcalc/REMD/AD/para/para_s_REMD_GBSA_ff99_2012-07-24.sh

dir=~/seminars/GS/2012-09-05

source ${parafile}

tiffheight=$( echo ${tiffwidth} \* 320 | bc )
tiffheight=$( echo ${tiffheight} / 390 | bc )

echo width=${tiffwidth}
echo height=${tiffheight}

paraRfil=${dir}/fig/para_fig_REMD_2Dpmf_GBSA_AD_2012-08-27.R

cat <<EOF > ${paraRfil}
pname<-"${pname}"

numEX<-"${numEX}"

TLbase<-"${TLbase}"

level <- seq(0,20,2)

ff<-"${ff}"

name.out <- paste("${dir}/tiff/fig_REMD_2Dpmf_GBSA_AD_2012-08-27",sep="")

tiffwidth <- ${tiffwidth}

tiffheight <- ${tiffheight}

source("${dir}/graph/graph_REMD_2Dpmf_GBSA_AD_2012-08-27.R")
EOF

Rscript ${paraRfil}; echo ${paraRfil}


