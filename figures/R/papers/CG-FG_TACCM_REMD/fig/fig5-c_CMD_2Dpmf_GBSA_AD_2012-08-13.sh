
dir=~/calspa/refcalc/CMD/AD

source ${dir}/para/para_s_CMD_T300_ff99SB_2012-07-24.sh

paraRfil=~/papers/CG-FG_TACCM_REMD/fig/para_fig5-b_CMD_2Dpmf_GBSA_AD_2012-08-13.R

cat <<EOF > ${paraRfil}
pname<-"${pname}"

T<-"${T}"

level <- seq(0,20,2)

ff<-"${ff}"

width <- "0.3"

title <- NULL

name.out <- paste("~/papers/CG-FG_TACCM_REMD/tiff/fig5-c_CMD_2Dpmf_GBSA_AD_2012-08-13",sep="")

source("~/papers/CG-FG_TACCM_REMD/graph/graph_CMD_2Dpmf_GBSA_AD_2012-08-13.R")

EOF

Rscript ${paraRfil}; echo ${paraRfil}
