#!~/bin/sh

opt=(dummy name parabasefile phsiflag )
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

height=20 
width1=1 
numuene=4 
n=20
nx=4
ny=5

dir=~/calspa/TACCM_CGAAREMD/AD
dirs=~/calspa/TACCM_CGAAREMD/AD

source ${parabasefile}

height=${height} 

AACGflag=CG

paraRfil=${dirs}/graph/para_fig_pmf1D_comp_vpara_e_CG-FG_NH_2012-11-11.sh

cat <<EOF > ${paraRfil}
nparameters <- ${nparameters}
tau <- NULL
mZ <- NULL
TZs <- NULL
pname1 <- NULL
numEX1 <- NULL
TLbase1 <- NULL
TAA <- NULL
TCG <- NULL
TZs <- NULL
KZAAo <- NULL
KZCGo <- NULL
AACG <- NULL

EOF

for  n in `seq 1 ${nparameters}`; do
    cat <<EOF >> ${paraRfil}
  TAA[$n]<-"${TAA_MuSTARMD[$n]}"
  TCG[$n]<-"${TCG_MuSTARMD[$n]}"
  TZs[$n] <- "${TZs_MuSTARMD[$n]}"
  tau[$n] <- "${tau_MuSTARMD[$n]}"
  mZ[$n]  <- "${mZ_MuSTARMD[$n]}"
#  nKZAA[$n]<-${nKZAA_MuSTARMD[$n]}
#  nKZCG[$n]<-${nKZCG_MuSTARMD[$n]}
#  numRE[$n]<-${numRE_MuSTARMD[$n]}
  pname1[$n]<-"${pname_MuSTARMD[$n]}"
  numEX1[$n]<-"${numEX_MuSTARMD[$n]}"
  TLbase1[$n]<-"${TLbase_MuSTARMD[$n]}"

  AACG[$n] <- "${AACGflag}"
  KZAAo[$n] <- "${KZAAo_MuSTARMD[$n]}"
  KZCGo[$n] <- "${KZCGo_MuSTARMD[$n]}"
EOF
done

cat <<EOF >> ${paraRfil}
width <- "${width1}"

numuene <- "${numuene}"

height <- ${height}

phsiflag <- "${phsiflag}"
phsi <- ${phsi}
name.title <- paste(AACG,sep="")

num <- ${n}
numx <- ${nx}
numy <- ${ny}
EOF

cat <<EOF >> ${paraRfil}

ff<-"${ffREMD}"

pname2<-"${pnameREMD}"

numEX2<-"${numEXREMD}"

TLbase2<-"${TLbaseREMD}"

name.out <- "${dirs}/tiff/fig_pmf1D_comp_vpara_e_CG-FG_NH_2012-11-11_${name}_${phsiflag}_supplement"

source("~/calspa/TACCM_CGAAREMD/AD/graph/graph_pmf1D_comp_vpara_e_CG-FG_NH_2012-11-11_supplement.R")
EOF
Rscript ${paraRfil}; echo ${paraRfil}

display ${dirs}/tiff/fig_pmf1D_comp_vpara_e_CG-FG_NH_2012-11-11_${name}_${phsiflag}_supplement.tiff &
