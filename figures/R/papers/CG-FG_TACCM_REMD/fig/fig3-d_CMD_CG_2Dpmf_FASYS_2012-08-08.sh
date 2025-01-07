
parafile=~/calspa/refcalc/CMD/FiveAtomSys/para/para_e_CMD_NH_2012-08-20_FG7_CG4_TL=1000000.sh

dir=~/calspa/refcalc/CMD/FiveAtomSys/

source ${parafile}

paraRfil=~/papers/CG-FG_TACCM_REMD/fig/para_fig3-c_CMD_CG_2Dpmf_FASYS_2012-08-08.R

i=2

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

name.out <- "~/papers/CG-FG_TACCM_REMD/tiff/fig3-d_CMD_CG_2Dpmf_FASYS_2012-08-08"

source("~/papers/CG-FG_TACCM_REMD/graph/graph_CMD_2Dpmf_FASYS_2012-08-08.R")

EOF

Rscript ${paraRfil}; echo ${paraRfil}

