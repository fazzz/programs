#!~/bin/sh

opt=(dummy phsiflag )
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

parafile1=~/calspa/TACCM_CGAAREMD/AD/para/2012-09-17/para_e_CG-FG_NR=4_wrefd0.1_TZ=350_fq=10ps_99SB_KZmax=1000.sh
parafile2=~/calspa/refcalc/REMD/AD/para/para_s_REMD_vac_ff99SB_2012-07-24.sh
height=20 
width1=1 
numuene=4 
n=20
nx=4
ny=5

dir=~/calspa/TACCM_CGAAREMD/AD
dirs=~/papers/CG-FG_TACCM_REMD

source ${parafile1}

height=${height} 

TZ=${TZs[1]}

AACGflag=CG

paraRfil=${dirs}/graph/para_1Dpmf_comp_MUSTERMD-TREMD_2012-09-25_${pname}_TZ=${TZ}_AA_2012-09-13.R

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

pname1<-"${pname}"

numEX1<-"${numEX}"

TLbase1<-"${TLbase}"

AACG <- "${AACGflag}"
width <- "${width1}"

KZAAo <- "${KZAAo[2]}"
KZCGo <- "${KZCGo[2]}"

numuene <- "${numuene}"

height <- ${height}

phsiflag <- "${phsiflag}"
phsi <- ${phsi}
name.title <- paste(AACG,sep="")

num <- ${n}
numx <- ${nx}
numy <- ${ny}
EOF

source ${parafile2}

cat <<EOF >> ${paraRfil}

ff<-"${ff}"

pname2<-"${pname}"

numEX2<-"${numEX}"

TLbase2<-"${TLbase}"

name.out <- "~/papers/CG-FG_TACCM_REMD/tiff/fig5_comp_MuSTARMD-TREMD_AD_1Dpmf_${phsiflag}_2012-10-03"

source("${dirs}/graph/graph_1Dpmf_comp_MUSTERMD-TREMD_2012-10-03.R")
EOF
Rscript ${paraRfil}; echo ${paraRfil}
