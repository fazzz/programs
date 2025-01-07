TAA_MuSTARMD<-"300"
TCG_MuSTARMD<-"300"

TZ_MuSTARMD <- "1000"

tau_MuSTARMD <- "1.0"

mZ_MuSTARMD <-"50000.00"  

pname_MuSTARMD <- "T1_FG7-CG4"

numEX_MuSTARMD<-"4000"

fq_MuSTARMD<-"10"

TLbase_MuSTARMD<-"10"

width_MuSTARMD <- "0.3"

T_CMD<-"300"

tau_CMD<-c(  "1.0" )

ffname_CMD_FG<-"FG7"
ffname_CMD_CG<-"CG4"

TL_CMD<-"1000000"

width_CMD <- "0.3"

level <- seq(0,10,1)

name.title <- NULL

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

file.name <- "~/papers/CG-FG_TACCM_REMD/tiff/fig2_a-b-c-d_comp_MuSTARMD-CMD_FASYS_2Dpmf_2012-10-03.tiff"
tiff(file.name,width=550,height=400)

par(oma = c(0,0,0,0) )
source("~/papers/CG-FG_TACCM_REMD/FillConMap.R")
source("~/papers/CG-FG_TACCM_REMD/FillConBar.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConBar.R")
						
nf <- layout(matrix(c(1,2,5,3,4,5),2,3,byrow=TRUE),c(3,3,1),c(1,1))

par(mar =c(5,5,2,2) )
par(cex.axis=2.0)
par(cex.lab=2.0)
label.x <- expression(paste(phi))
label.y <- expression(paste(psi))

title=name.title
	
dir0 <- "/home/yamamori/calspa/TACCM_CGAAREMD/FiveAtomSysf"
dirbase <- paste(dir0,"/e_CG-FG_NH_2012-08-14_2",sep="")
AACG <- "CG"
KZAAo <- "0"
KZCGo <- "1000"
name <- paste(dirbase,"/tau=",tau_MuSTARMD,"/mZ=",mZ_MuSTARMD,"/TZ=",TZ_MuSTARMD,"/",pname_MuSTARMD,"/nEX=",numEX_MuSTARMD,"/fq=",fq_MuSTARMD,"/pmf/pmf_TAA=",TAA_MuSTARMD,"_TCG_",TCG_MuSTARMD,"_TZ_",TZ_MuSTARMD,"_",width_MuSTARMD,"_",KZAAo,"_",KZCGo,"_",AACG,"_4",sep="")
felwFillConMap(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

AACG <- "AA"
KZAAo <- "1000"
KZCGo <- "0"
name <- paste(dirbase,"/tau=",tau_MuSTARMD,"/mZ=",mZ_MuSTARMD,"/TZ=",TZ_MuSTARMD,"/",pname_MuSTARMD,"/nEX=",numEX_MuSTARMD,"/fq=",fq_MuSTARMD,"/pmf/pmf_TAA=",TAA_MuSTARMD,"_TCG_",TCG_MuSTARMD,"_TZ_",TZ_MuSTARMD,"_",width_MuSTARMD,"_",KZAAo,"_",KZCGo,"_",AACG,"_4",sep="")
felwFillConMap(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

dir0 <- "/home/yamamori/calspa/refcalc/CMD/FiveAtomSys"
dirbase <- paste(dir0,"/e_CMD_NH_2012-08-20",sep="")
name <- paste(dirbase,"/ff=",ffname_CMD_FG,"/tau=",tau_CMD,"/TL=",TL_CMD,"/anl/pmf_FASYSv_T",T_CMD,"_",width_CMD,sep="")
felwFillConMap(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

dir0 <- "/home/yamamori/calspa/refcalc/CMD/FiveAtomSys"
dirbase <- paste(dir0,"/e_CMD_NH_2012-08-20",sep="")
name <- paste(dirbase,"/ff=",ffname_CMD_CG,"/tau=",tau_CMD,"/TL=",TL_CMD,"/anl/pmf_FASYSv_T",T_CMD,"_",width_CMD,sep="")
felwFillConMap(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

felwFillConBar(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

