
TAA<-"300"
TCG<-"300"

TZs <- c( "500" )
nTZs <- length(TZs)

tau<-c(  "1.0" )
ntau<-length(tau)

mZ<-c(  "50000.00" )
nmZ<-length(mZ)

KZAA<-c( "1000", "1000", "1000", "1000", "500", "250", "125", "100"  ) 
KZCG<-c( "100",  "125",  "250",  "500", "1000", "1000", "1000",  "1000"  )

nKZAA<-length(KZAA)
nKZCG<-length(KZCG)
numRE<-length(KZAA)

pname<-"T1"

numEX<-"10000"

AACG <- "AA"
width <- "0.1"

KZAAo <- "100" 
KZCGo <- "0"   

#####################
## KZAAo <- "1000" ##
## KZCGo <- "0"    ##
#####################

level <- seq(0,25,1)

name.title <- NULL

name.out <- "~/summary/CG-FG_TACCM_REMD/R/fig3-a_TACCM_REMD_FG_2Dpmf_AD_2012-05-24"

source("~/summary/CG-FG_TACCM_REMD/R/graph_TACCM_REMD_FG_2Dpmf_AD_2012-05-24.R")

#OutTiff(name.out,w_val=390,h_val=320)
