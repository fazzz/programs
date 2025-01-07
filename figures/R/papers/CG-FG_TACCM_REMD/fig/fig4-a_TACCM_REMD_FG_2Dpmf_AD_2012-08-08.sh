
dir=~/papers/CG-FG_TACCM_REMD/

parafile=~/calspa/TACCM_CGAAREMD/AD/para/para_e_CG-FG_NH_2012-07-27_wrefd0.1_TZ=750_fq=10ps_99SB.sh 

source ${parafile}

paraRfil=${dir}/fig/para_fig4-a_TACCM_REMD_FG_2Dpmf_AD_2012-08-08.R

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

AACG <- "AA"
width <- "0.3"

KZAAo <- "1000" 
KZCGo <- "0"   

level <- seq(0,20,2)

name.title <- NULL

name.out <- "~/papers/CG-FG_TACCM_REMD/tiff/fig4-a_TACCM_REMD_FG_2Dpmf_AD_2012-08-08"

source("~/papers/CG-FG_TACCM_REMD/graph/graph_TACCM_REMD_2Dpmf_AD_2012-08-08.R")

EOF

Rscript ${paraRfil}; echo ${paraRfil}
