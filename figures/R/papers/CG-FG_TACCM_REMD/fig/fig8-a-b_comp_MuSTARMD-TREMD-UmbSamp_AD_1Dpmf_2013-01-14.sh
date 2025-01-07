#!~/bin/sh

phsiflag=phi

parafile1=~/calspa/TACCM_CGAAREMD/AD/para/2012-11-02/para_e_CG-FG_NR=4_TZ=700_fq=10ps_99SB_KZmax=5000_weljd0.001_mZ=100.sh

parafile2=~/calspa/refcalc/REMD/AD/para/para_s_REMD_vac_ff99SB_2012-07-24.sh
parafileUmbSam=~/calspa/refcalc/UmbSam/AD/para/para_s_vac_ff99SB_K=10_2012-11-12.sh
height=20 
width1=1 
numuene=4 
n=20
nx=4
ny=5

dir=~/calspa/TACCM_CGAAREMD/AD
dirOUT=~/papers/CG-FG_TACCM_REMD

source ${parafile1}

height=${height} 

TZ=${TZs[1]}

AACGflag=CG

paraRfil=${dirOUT}/graph/para_comp_MuSTARMD-TREMD-UmbSamp_AD_1Dpmf_2013-01-14.R

cat <<EOF > ${paraRfil}
TAA<-"${TAA}"
TCG<-"${TCG}"

TZs <- c( "${TZ}" )
nTZs <- length(TZs)

tau<-c(  "1.0" )
ntau<-length(tau)

mZ<-c(  "${mZ[1]}" )
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
EOF

source ${parafileUmbSam}
TLns=`expr ${TLbase} / 1000`

cat <<EOF >> ${paraRfil}
ff2<-"${ff}"

pname3<-"${pname}"

TLns<-"${TLns}"

name.out <- "${dirOUT}/tiff/fig8-a_comp_MuSTARMD-TREMD-UmbSamp_AD_1Dpmf_2013-01-14"

source("${dirOUT}/graph/graph_1Dpmf_MuSTARMD_comp_TREMD_2DFEL_2013-01-14.R")
EOF

Rscript ${paraRfil}; echo ${paraRfil}

####################################################################################################

phsiflag=psi

parafile1=~/calspa/TACCM_CGAAREMD/AD/para/2012-11-02/para_e_CG-FG_NR=4_TZ=700_fq=10ps_99SB_KZmax=5000_weljd0.001_mZ=100.sh

parafile2=~/calspa/refcalc/REMD/AD/para/para_s_REMD_vac_ff99SB_2012-07-24.sh
parafileUmbSam=~/calspa/refcalc/UmbSam/AD/para/para_s_vac_ff99SB_K=10_2012-11-12.sh
height=20 
width1=1 
numuene=4 
n=20
nx=4
ny=5

dir=~/calspa/TACCM_CGAAREMD/AD
dirOUT=~/papers/CG-FG_TACCM_REMD

source ${parafile1}

height=${height} 

TZ=${TZs[1]}

AACGflag=CG

paraRfil=${dirOUT}/graph/para_comp_MuSTARMD-TREMD-UmbSamp_AD_1Dpmf_2013-01-14.R

cat <<EOF > ${paraRfil}
TAA<-"${TAA}"
TCG<-"${TCG}"

TZs <- c( "${TZ}" )
nTZs <- length(TZs)

tau<-c(  "1.0" )
ntau<-length(tau)

mZ<-c(  "${mZ[1]}" )
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
EOF

source ${parafileUmbSam}
TLns=`expr ${TLbase} / 1000`

cat <<EOF >> ${paraRfil}
ff2<-"${ff}"

pname3<-"${pname}"

TLns<-"${TLns}"

name.out <- "${dirOUT}/tiff/fig8-b_comp_MuSTARMD-TREMD-UmbSamp_AD_1Dpmf_2013-01-14"

source("${dirOUT}/graph/graph_1Dpmf_MuSTARMD_comp_TREMD_2DFEL_2013-01-14.R")
EOF

Rscript ${paraRfil}; echo ${paraRfil}
