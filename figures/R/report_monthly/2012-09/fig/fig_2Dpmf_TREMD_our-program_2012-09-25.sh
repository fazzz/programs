#!~/bin/sh

opt=(dummy size  )
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

###################################################
# width_size=${size}				  #
# height_size=`echo  ${width_size} / 390 | bc` 	  #
# height_size=`echo  ${height_size} \* 320 | bc`  #
###################################################

width_size=390
height_size=320

height=10
dwidth=1
parafile=~/calspa/refcalc/REMD/AD/para/para_s_REMD_vac_ff99SB_2012-07-24.sh

dir=~/calspa/refcalc/REMD/AD
direqubase=~/calspa/refcalc/REMD/AD/s_REVAC_2012-07-23_${ff}/

source ${parafile}

paraRfil=~/Report/2012-09/graph/para_2Dpmf_TREMD_our-program_2012-09-25.R

cat <<EOF > ${paraRfil}
pname<-"${pname}"

numEX<-"${numEX}"

TLbase<-"${TLbase}"

level <- seq(0,${height},${dwidth})

ff<-"${ff}"

name.out <- paste("~/Report/2012-09/tiff/pmf_2D_TREMD_our-program_2012-09-25",sep="")

size <- ${size}

source("~/Report/2012-09/graph/graph_2Dpmf_TREMD_our-program_2012-09-25.R")
EOF
Rscript ${paraRfil}; echo ${paraRfil}
