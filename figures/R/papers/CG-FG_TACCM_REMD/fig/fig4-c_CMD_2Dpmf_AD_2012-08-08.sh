
dir=~/calspa/refcalc/CMD/AD

source ${dir}/para/para_s_CMD_T300_ff99SB_2012-07-24.sh

paraRfil=~/papers/CG-FG_TACCM_REMD/fig/para_fig4-c_CMD_2Dpmf_AD_2012-08-08.R

cat <<EOF > ${paraRfil}
pname<-"${pname}"

T<-"${T}"

level <- seq(0,20,2)

ff<-"${ff}"

width <- "0.2"

title <- NULL

name.out <- paste("~/papers/CG-FG_TACCM_REMD/tiff/fig4-c_CMD_2Dpmf_AD_2012-08-08",sep="")

source("~/papers/CG-FG_TACCM_REMD/graph/graph_CMD_2Dpmf_AD_2012-08-08.R")

EOF

Rscript ${paraRfil}; echo ${paraRfil}
