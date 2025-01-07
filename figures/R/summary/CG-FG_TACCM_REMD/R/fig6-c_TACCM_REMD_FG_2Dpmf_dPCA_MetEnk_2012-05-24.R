
TAA<-"300"
TCG1<-"370"
TCG2<-"370"

TZs <- c( "1000" )
nTZs <- length(TZs)

tau<-c(  "1.0" )
ntau<-length(tau)

ep<-"0.2"
nb<-"3"
cutoff<-"4.7"

mZ<-c(  "50000.00" )
nmZ<-length(mZ)

nKZAA<-8
nKZCG1<-8
nKZCG2<-8
numRE<-8

pname<-"2CG11"

numEX<-"10000"

AACG <- "CG2"
width <- "0.01"
#width <- "0.005"

KZAAo <- "0"
KZCG1o <- "0"
KZCG2o <- "100"

#level <- seq(-60,0,6)
level <- seq(-200,0,10)

state <- "1FG2CG"

name.title <- NULL

name.out <- "~/summary/CG-FG_TACCM_REMD/R/fig6-a_TACCM_REMD_FG_2Dpmf_dPCA_MetEnk_2012-05-24"

source("~/summary/CG-FG_TACCM_REMD/R/graph_TACCM_REMD_2Dpmf_dPCA_MetEnk_2012-05-24.R")

OutTiff(name.out,w_val=390,h_val=320)
