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

dir=~/calspa/refcalc/CMD/AD

source ${dir}/para/para_s_CMD_T300_ff99SB_2012-07-24.sh

tiffheight=$( echo ${tiffwidth} \* 320 | bc )
tiffheight=$( echo ${tiffheight} / 390 | bc )

echo width=${tiffwidth}
echo height=${tiffheight}

paraRfil=~/seminars/GS/2012-09-05/fig/para_fig4-c_CMD_2Dpmf_AD_2012-08-27.R

cat <<EOF > ${paraRfil}
pname<-"${pname}"

T<-"${T}"

level <- seq(0,20,2)

ff<-"${ff}"

width <- "0.2"

title <- NULL

name.out <- paste("~/seminars/GS/2012-09-05/tiff/fig_CMD_2Dpmf_AD_2012-08-27",sep="")

tiffwidth <- ${tiffwidth}

tiffheight <- ${tiffheight}

source("~/seminars/GS/2012-09-05/graph/graph_CMD_2Dpmf_AD_2012-08-27.R")

EOF

Rscript ${paraRfil}; echo ${paraRfil}
