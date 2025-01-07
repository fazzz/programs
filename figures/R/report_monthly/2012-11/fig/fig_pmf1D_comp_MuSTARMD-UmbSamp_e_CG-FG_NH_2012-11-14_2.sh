#!~/bin/sh

opt=(dummy tiffheight tiffwidth )
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

phsiflag=psi
num1=4

height=20 
width1=1 
numuene=4 
n=3
nx=4
ny=3

dir=~/calspa/TACCM_CGAAREMD/AD
dirout=~/Report/2012-11/

parafileUmbSam=~/calspa/refcalc/UmbSam/AD/para/para_s_vac_ff99SB_K=10_2012-11-12.sh
source ${parafileUmbSam}

ffUmbSam=${ff}
pnameUmbSam=${pname}
TLns=`expr ${TLbase} / 1000`

parafileMuSTARMD=( dummy  ~/calspa/TACCM_CGAAREMD/AD/para/2012-11-02/para_e_CG-FG_NR=4_TZ=500_fq=10ps_99SB_KZmax=1000_weljd0.001_mZ=100.sh ~/calspa/TACCM_CGAAREMD/AD/para/2012-11-02/para_e_CG-FG_NR=4_TZ=500_fq=10ps_99SB_KZmax=2500_weljd0.001_mZ=100.sh ~/calspa/TACCM_CGAAREMD/AD/para/2012-11-02/para_e_CG-FG_NR=4_TZ=500_fq=10ps_99SB_KZmax=5000_weljd0.001_mZ=100.sh )

nparameters=`expr ${#parafileMuSTARMD[*]} - 1`

tau_MuSTARMD=( dummy )
mZ_MuSTARMD=( dummy )
TZs_MuSTARMD=( dummy )
pname_MuSTARMD=( dummy )
numEX_MuSTARMD=( dummy )
TLbase_MuSTARMD=( dummy )
TAA_MuSTARMD=( dummy )
TCG_MuSTARMD=( dummy )
TZs_MuSTARMD=( dummy )
KZAAo_MuSTARMD=( dummy )
KZCGo_MuSTARMD=( dummy )

for n in `seq 1 ${nparameters}`; do
    source ${parafileMuSTARMD[$n]}

    tau_MuSTARMD[$n]=${tau[1]}
    mZ_MuSTARMD[$n]=${mZ[1]}
    TZs_MuSTARMD[$n]=${TZs[1]}
    pname_MuSTARMD[$n]=${pname}
    numEX_MuSTARMD[$n]=${numEX}
    TLbase_MuSTARMD[$n]=${TLbase}
    TAA_MuSTARMD[$n]=${TAA}
    TCG_MuSTARMD[$n]=${TCG}
    KZAAo_MuSTARMD[$n]=${KZAAo[2]}
    KZCGo_MuSTARMD[$n]=${KZCGo[2]}
done

TZs_MuSTARMD=( dummy 350 500 600 700 )
nTz=`expr ${#TZs_MuSTARMD[*]} - 1`

mZ_MuSTARMD=( dummy "100.00" "1000.00" "10000.00" )
nmZ=`expr ${#mZ_MuSTARMD[*]} - 1`

height=${height} 

AACGflag=CG

paraRfil=${dirout}/graph/para_fig_pmf1D_comp_MuSTARMD-UmbSamp_Kzdependence_e_CG-FG_NH_2012-11-14.sh

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
  tau[$n] <- "${tau_MuSTARMD[$n]}"
  mZ[$n]  <- "${mZ_MuSTARMD[$n]}"
  pname1[$n]<-"${pname_MuSTARMD[$n]}"
  numEX1[$n]<-"${numEX_MuSTARMD[$n]}"
  TLbase1[$n]<-"${TLbase_MuSTARMD[$n]}"

  AACG[$n] <- "${AACGflag}"
  KZAAo[$n] <- "${KZAAo_MuSTARMD[$n]}"
  KZCGo[$n] <- "${KZCGo_MuSTARMD[$n]}"
EOF
done

for  n in `seq 1 ${nTz}`; do
    cat <<EOF >> ${paraRfil}
  TZs[$n] <- "${TZs_MuSTARMD[$n]}"
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

ff<-"${ffUmbSam}"

pname2<-"${pnameUmbSam}"

TLns<-"${TLns}"

angle <- ${num1}

N <- 3

name.out <- "${dirout}/tiff/fig_pmf1D_comp_MuSTARMD-UmbSamp_e_CG-FG_NH_2012-11-14_${phsiflag}_${num1}"

tiffwidth=${tiffwidth}
tiffheight=${tiffheight}

source("~/Report/2012-11/graph/graph_pmf1D_comp_MuSTARMD-UmbSamp_vpara_e_CG-FG_NH_2012-11-14.R")
EOF

Rscript ${paraRfil}; echo ${paraRfil}

