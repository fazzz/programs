nparameters <- 3
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

  TAA[1]<-"300"
  TCG[1]<-"300"
  tau[1] <- "1.0"
  mZ[1]  <- "100.00"
  pname1[1]<-"SB_KZMAX=1000_NR4_woeljd0.001"
  numEX1[1]<-"1000"
  TLbase1[1]<-"10"

  AACG[1] <- "CG"
  KZAAo[1] <- "0"
  KZCGo[1] <- "1000"
  TAA[2]<-"300"
  TCG[2]<-"300"
  tau[2] <- "1.0"
  mZ[2]  <- "1000.00"
  pname1[2]<-"SB_KZMAX=2500_NR4_woeljd0.001"
  numEX1[2]<-"1000"
  TLbase1[2]<-"10"

  AACG[2] <- "CG"
  KZAAo[2] <- "0"
  KZCGo[2] <- "2500"
  TAA[3]<-"300"
  TCG[3]<-"300"
  tau[3] <- "1.0"
  mZ[3]  <- "10000.00"
  pname1[3]<-"SB_KZMAX=5000_NR4_woeljd0.001"
  numEX1[3]<-"1000"
  TLbase1[3]<-"10"

  AACG[3] <- "CG"
  KZAAo[3] <- "0"
  KZCGo[3] <- "5000"
  TZs[1] <- "350"
  TZs[2] <- "500"
  TZs[3] <- "600"
  TZs[4] <- "700"
width <- "1"

numuene <- "4"

height <- 20

phsiflag <- "psi"
phsi <- 
name.title <- paste(AACG,sep="")

num <- 4
numx <- 4
numy <- 3

ff<-"ff99SB"

pname2<-"Umb_Nbin=12x12_K=10"

TLns<-"10"

angle <- 4

N <- 3

name.out <- "/home/yamamori/Report/2012-11//tiff/fig_pmf1D_comp_MuSTARMD-UmbSamp_e_CG-FG_NH_2012-11-14_psi_4"

tiffwidth=650
tiffheight=340

source("~/Report/2012-11/graph/graph_pmf1D_comp_MuSTARMD-UmbSamp_vpara_e_CG-FG_NH_2012-11-14.R")
