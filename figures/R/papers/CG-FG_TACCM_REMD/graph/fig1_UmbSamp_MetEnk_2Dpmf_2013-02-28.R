TAA<-"300"
TCG1<-"370"
TCG2<-"370"

TZs <-"1400"
tau<-"1.0"

ep<-"0.2"
nb<-"3"
cutoff<-"4.7"

mZ<-"50000.00"

nKZAA<-8
nKZCG1<-8
nKZCG2<-8
numRE<-8

pname<-"2CG21_KAA=1000"

freq<-"1"

AACG <- "AA"
width <- "0.01"

KZAAo <- "1000"
KZCG1o <- "0"
KZCG2o <- "0"

level <- seq(0,5,0.5)

state <- "1FG2CG"

TS <- "10000"

TLbase <- "1"

MyColor <- function(n,alpha=1)
{			
    if ((n <- as.integer(n[1L])) > 0) {
        j <- n%/%3		
        k <- n%/%3		
        i <- n - j - k		
        c(if (i > 0) hsv(h = seq.int(from = 40/60, to = 25/60,
            length.out = i), alpha = alpha), if (j > 0) hsv(h = seq.int(from = 23/60,
            to = 11/60, length.out = j), alpha = alpha), if (k > 
            0) hsv(h = seq.int(from = 8/60, to = 0/60, length.out = k-1),
            alpha = alpha, s = seq.int(from = 1, to = 0.9, length.out = k-1),
            v = 1),hsv(0,0,1))
    }			
    else character(0L)	
}			

file.name <- "~/papers/CG-FG_TACCM_REMD/eps/fig1_UmbSamp_MetEnk_2Dpmf_2013-02-28.eps"
#tiff(file.name,width=800,height=600)
postscript(file.name,width=3.2,height=2.4,horizontal=FALSE,onefile=FALSE,paper="special")

xrange <- c(0,0.3,6)
yrange <- c(0,0.3,6)

par(oma = c(0,0,0,0) )
source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange.R")
source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange_wy.R")
source("~/papers/CG-FG_TACCM_REMD/FillConBar.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap_wrange.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap_wrange_wy.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConBar.R")

nf <- layout(matrix(c(1,2),1,2,byrow=TRUE),c(3.0,1),c(1))

par(cex.axis=1.0)
par(cex.lab=1.0)
label.x <- "d1"
label.y <- "d2"

name.title <- NULL

title=name.title
	
dir0 <- "/home/yamamori/calspa/TACCM_CGAAREMD/MetEnk/"
dirbase <- paste(dir0,"/e_CG-FG_CMD_NH_2012-06-04/",state,sep="")

nTZs<-1

name <- paste(dirbase,"/tau=",tau,"/mZ=",mZ,"/ep=",ep,"/cutoff=",cutoff,"/TZ=",TZs,"/",pname,"/freq=",freq,"/pmf/pmf_TAA=",TAA,"_TCG1_",TCG1,"_TCG2_",TCG2,"_TZ_",TZs,"_",width,"_",KZAAo,"_",KZCG1o,"_",KZCG2o,"_",AACG,"_bo10000ps_2",sep="")

cat(name)

felwFillConMapwrangewy(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,mar=c(4,4,1,1))

par(mar=c(4,2,1,1))

felwFillConBar(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

