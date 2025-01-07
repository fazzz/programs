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

parafile=~/calspa/refcalc/CMD/FiveAtomSys/para/para_e_CMD_NH_2012-08-20_FG7_CG4_TL=1000000.sh

dir=~/calspa/refcalc/CMD/FiveAtomSys/

source ${parafile}

tiffheight=$( echo ${tiffwidth} \* 320 | bc )
tiffheight=$( echo ${tiffheight} / 390 | bc )

echo width=${tiffwidth}
echo height=${tiffheight}

paraRfil=~/papers/CG-FG_TACCM_REMD/fig/para_fig_CMD_FG_2Dpmf_FASYS_2012-08-27.R

i=1

cat <<EOF > ${paraRfil}
T<-"${T}"

tau<-c(  "1.0" )
ntau<-length(tau)

ffname<-"${ffname[$i]}"
nffname<-1

TL<-"${TLbase}"

width <- "0.3"

level <- seq(0,10,1)

name.title <- NULL

name.out <- "~/seminars/GS/2012-09-05/tiff/fig_CMD_FG_2Dpmf_FASYS_2012-08-27"

tiffwidth <- ${tiffwidth}

tiffheight <- ${tiffheight}

source("~/seminars/GS/2012-09-05/graph/graph_CMD_2Dpmf_FASYS_2012-08-27.R")

EOF

Rscript ${paraRfil}; echo ${paraRfil}

