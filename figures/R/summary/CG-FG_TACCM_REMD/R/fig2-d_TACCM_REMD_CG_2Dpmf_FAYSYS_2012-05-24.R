
TAA<-"300"
TCG<-"370"

TZs <- c( "1000" )
nTZs <- length(TZs)

tau<-c(  "1.0" )
ntau<-length(tau)

mZ<-c(  "50000.00" )
nmZ<-length(mZ)

KZAA<-c( "1000", "1000", "1000", "1000", "1000", "500",   "250",  "125"  ) 
KZCG<-c(  "100",  "125",  "250",  "500", "1000", "1000", "1000", "1000"  )

nKZAA<-length(KZAA)
nKZCG<-length(KZCG)
numRE<-length(KZAA)

pname<-"T6"

numEX<-"1000"

AACG <- "AA"
width <- "0.1"

KZAAo <- "3000"
KZCGo <- "0"

level <- seq(0,20,1)

name.title <- NULL

#source("~/calspa/TACCM_CGAAREMD/FiveAtomSysd/R/graph_pmf2D_TACCM_REMD_e_CG-FG_NH_2012-05-01.R")
name.out <- "~/summary/CG-FG_TACCM_REMD/R/fig2-d_TACCM_REMD_FG_2Dpmf_FAYSYS_2012-05-24"
source("~/summary/CG-FG_TACCM_REMD/R/graph_TACCM_REMD_2Dpmf_FAYSYS_2012-05-24.R")

#OutTiff(name.out,w_val=390,h_val=300)
