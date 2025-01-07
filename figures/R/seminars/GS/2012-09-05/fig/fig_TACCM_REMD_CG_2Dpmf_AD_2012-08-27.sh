#!~/bin/sh

opt=(dummy name tiffwidth parafile )
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

#parafile=~/calspa/TACCM_CGAAREMD/AD/para/para_e_CG-FG_NH_2012-07-27_wrefd0.1_TZ=750_fq=10ps_99SB.sh 

source ${parafile}

tiffheight=$( echo ${tiffwidth} \* 320 | bc )
tiffheight=$( echo ${tiffheight} / 390 | bc )

echo width=${tiffwidth}
echo height=${tiffheight}

paraRfil=${dir}/fig/para_fig_TACCM_REMD_CG_2Dpmf_AD_2012-08-27.R

cat <<EOF > ${paraRfil}

TAA<-"${TAA}"
TCG<-"${TCG}"

TZ <-"${TZs[1]}"

tau<-"${tau[1]}"

mZ <-"${mZ[1]}"  

pname<-"${pname}"

numEX<-"${numEX}"

fq<-"${TLbase}"

TLbase<-"${TLbase}"

AACG <- "CG"
width <- "0.3"

KZAAo <- "0" 
KZCGo <- "1000"   

level <- seq(0,20,2)

name.title <- NULL

name.out <- "${dir}/tiff/fig_TACCM_REMD_CG_2Dpmf_AD_2012-08-27_${name}"

tiffwidth <- ${tiffwidth}

tiffheight <- ${tiffheight}

source("${dir}/graph/graph_TACCM_REMD_2Dpmf_AD_2012-08-27.R")

EOF

Rscript ${paraRfil}; echo ${paraRfil}
