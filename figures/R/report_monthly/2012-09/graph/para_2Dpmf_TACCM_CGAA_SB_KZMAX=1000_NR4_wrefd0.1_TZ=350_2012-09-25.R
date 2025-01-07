TAA<-"300"
TCG<-"300"

TZs <- c( "350" )
nTZs <- length(TZs)

tau<-c(  "1.0" )
ntau<-length(tau)

mZ<-c(  "50000.00" )
nmZ<-length(mZ)

nKZAA<-4
nKZCG<-4
numRE<-4

pname<-"SB_KZMAX=1000_NR4_wrefd0.1"

numEX<-"1000"

TLbase<-"10"

AACG <- "CG"
width <- "0.3"

KZAAo <- "0"
KZCGo <- "1000"

numuene <- "4"

level <- seq(0,10,1)

name.title <- paste(AACG,sep="")

name.out <- "~/Report/2012-09/tiff/pmf_2D_TACCM_CGAA_our-program_SB_KZMAX=1000_NR4_wrefd0.1_TZ=350_2012-09-25"

size <- 1.0

source("~/Report/2012-09/graph/graph_2Dpmf_TACCM_CGAA_our-program_2012-09-25.R")
