fact.x=1
fact.y=1

level <- seq(0,20,2)

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

file.name <- "~/Report/2013-08/eps/fig_a-h_40ns-TAMD_2Dpmf_AD_2013-08-17.eps"

postscript(file.name,width=6.0,height=9.3,horizontal=FALSE,onefile=FALSE,paper="special")

xrange <- c(-1.0,1.0,4)
xrange.axis <- c(-1.0,1.0,4)
yrange <- c(-1.0,1.0,4)
yrange.axis <- c(-1.0,1.0,4)

par(mar = c(0.0,0.0,0.0,0.0) )
par(oma = c(0,0,0,0) )

source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange_worwoaxis.R")
source("~/papers/CG-FG_TACCM_REMD/FillConBar.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap_wrange_worwoaxis.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConBar_wmai.R")

source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange_worwoaxis.R")
source("~/papers/CG-FG_TACCM_REMD/FillConBar.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap_wrange_worwoaxis_wpoints.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConBar_wmai.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMapbox_wpoints_2.R")

nf <- layout(matrix(c(1,2,9,3,4,9,5,6,9,7,8,9),4,3,byrow=TRUE),c(2.9,2.1,1.1),c(2.5,2.0,2.0,2.8))

par(cex.axis=1.0)
par(cex.lab=1.0)

label.x <- ""
label.y <- ""

title=name.title
T_TAMD<-"300"
TB_TAMD<-"750"
tau_TAMD<-"1.0"
KZ_TAMD<-"1000"
mZ_TAMD<-"100.00"
width_TAMD <- "0.3"
ff_TAMD <- "99SB"

dir0 <- "/home/yamamori/calspa/TACCM/AD"
dirbase <- paste(dir0,"/e_TACCM_NH_2012-07-31_",ff_TAMD,sep="")
name <- paste(dirbase,"/tau=",tau_TAMD,"/TB=",TB_TAMD,"/KZ=",KZ_TAMD,"/mZ=",mZ_TAMD,"/pmf/pmf_T=",T_TAMD,".Zhist_0.2",sep="")

cat(name)
	
felwFillConMapwrangeworwoaxis(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,xaflag="f",mai=c(0.0,0.8,0.5,0.0),norm="T")
T_TAMD<-"300"
TB_TAMD<-"750"
tau_TAMD<-"1.0"
KZ_TAMD<-"1000"
mZ_TAMD<-"10000.00"
width_TAMD <- "0.3"
ff_TAMD <- "99SB"

dir0 <- "/home/yamamori/calspa/TACCM/AD"
dirbase <- paste(dir0,"/e_TACCM_NH_2012-07-31_",ff_TAMD,sep="")
name <- paste(dirbase,"/tau=",tau_TAMD,"/TB=",TB_TAMD,"/KZ=",KZ_TAMD,"/mZ=",mZ_TAMD,"/pmf/pmf_T=",T_TAMD,".Zhist_0.2",sep="")

cat(name)

yrange <- c(-1.0,0.5,3)
yrange.axis <- c(-1.0,0.5,3)
	
felwFillConMapwrangeworwoaxis(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,xaflag="f",yaflag="f",mai=c(0.0,0.0,0.5,0.1),norm="T")
T_TAMD<-"300"
TB_TAMD<-"750"
tau_TAMD<-"1.0"
KZ_TAMD<-"1000"
mZ_TAMD<-"100.00"
width_TAMD <- "0.3"
ff_TAMD <- "99SB"

dir0 <- "/home/yamamori/calspa/TACCM/AD"
dirbase <- paste(dir0,"/e_TACCM_NH_2012-07-31_",ff_TAMD,sep="")
name <- paste(dirbase,"/tau=",tau_TAMD,"/TB=",TB_TAMD,"/KZ=",KZ_TAMD,"/mZ=",mZ_TAMD,"/pmf/pmf_T=",T_TAMD,".Zhist_0.2",sep="")

cat(name)

yrange <- c(-1.0,0.5,3)
yrange.axis <- c(-1.0,0.5,3)
	
felwFillConMapwrangeworwoaxis(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,xaflag="f",mai=c(0.0,0.8,0.0,0.0),norm="T")
T_TAMD<-"300"
TB_TAMD<-"800"
tau_TAMD<-"1.0"
KZ_TAMD<-"1000"
mZ_TAMD<-"10000.00"
width_TAMD <- "0.3"
ff_TAMD <- "99SB"

dir0 <- "/home/yamamori/calspa/TACCM/AD"
dirbase <- paste(dir0,"/e_TACCM_NH_2012-07-31_",ff_TAMD,sep="")
name <- paste(dirbase,"/tau=",tau_TAMD,"/TB=",TB_TAMD,"/KZ=",KZ_TAMD,"/mZ=",mZ_TAMD,"/pmf/pmf_T=",T_TAMD,".Zhist_0.2",sep="")

cat(name)

yrange <- c(-1.0,0.5,3)
yrange.axis <- c(-1.0,0.5,3)
	
felwFillConMapwrangeworwoaxis(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,xaflag="f",yaflag="f",mai=c(0.0,0.0,0.0,0.1),norm="T")
T_TAMD<-"300"
TB_TAMD<-"1000"
tau_TAMD<-"1.0"
KZ_TAMD<-"1000"
mZ_TAMD<-"100.00"
width_TAMD <- "0.3"
ff_TAMD <- "99SB"

dir0 <- "/home/yamamori/calspa/TACCM/AD"
dirbase <- paste(dir0,"/e_TACCM_NH_2012-07-31_",ff_TAMD,sep="")
name <- paste(dirbase,"/tau=",tau_TAMD,"/TB=",TB_TAMD,"/KZ=",KZ_TAMD,"/mZ=",mZ_TAMD,"/pmf/pmf_T=",T_TAMD,".Zhist_0.2",sep="")

cat(name)

yrange <- c(-1.0,0.5,3)
yrange.axis <- c(-1.0,0.5,3)
	
felwFillConMapwrangeworwoaxis(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,xaflag="f",mai=c(0.0,0.8,0.0,0.0),norm="T")
T_TAMD<-"300"
TB_TAMD<-"1000"
tau_TAMD<-"1.0"
KZ_TAMD<-"1000"
mZ_TAMD<-"10000.00"
width_TAMD <- "0.3"
ff_TAMD <- "99SB"

dir0 <- "/home/yamamori/calspa/TACCM/AD"
dirbase <- paste(dir0,"/e_TACCM_NH_2012-07-31_",ff_TAMD,sep="")
name <- paste(dirbase,"/tau=",tau_TAMD,"/TB=",TB_TAMD,"/KZ=",KZ_TAMD,"/mZ=",mZ_TAMD,"/pmf/pmf_T=",T_TAMD,".Zhist_0.2",sep="")

cat(name)

yrange <- c(-1.0,0.5,3)
yrange.axis <- c(-1.0,0.5,3)
	
felwFillConMapwrangeworwoaxis(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,xaflag="f",yaflag="f",mai=c(0.0,0.0,0.0,0.1),norm="T")
T_TAMD<-"300"
TB_TAMD<-"1200"
tau_TAMD<-"1.0"
KZ_TAMD<-"1000"
mZ_TAMD<-"100.00"
width_TAMD <- "0.3"
ff_TAMD <- "99SB"

dir0 <- "/home/yamamori/calspa/TACCM/AD"
dirbase <- paste(dir0,"/e_TACCM_NH_2012-07-31_",ff_TAMD,sep="")
name <- paste(dirbase,"/tau=",tau_TAMD,"/TB=",TB_TAMD,"/KZ=",KZ_TAMD,"/mZ=",mZ_TAMD,"/pmf/pmf_T=",T_TAMD,".Zhist_0.2",sep="")

cat(name)

yrange <- c(-1.0,0.5,3)
yrange.axis <- c(-1.0,0.5,3)
	
felwFillConMapwrangeworwoaxis(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,mai=c(0.8,0.8,0.0,0.0),norm="T")
T_TAMD<-"300"
TB_TAMD<-"1200"
tau_TAMD<-"1.0"
KZ_TAMD<-"1000"
mZ_TAMD<-"10000.00"
width_TAMD <- "0.3"
ff_TAMD <- "99SB"

dir0 <- "/home/yamamori/calspa/TACCM/AD"
dirbase <- paste(dir0,"/e_TACCM_NH_2012-07-31_",ff_TAMD,sep="")
name <- paste(dirbase,"/tau=",tau_TAMD,"/TB=",TB_TAMD,"/KZ=",KZ_TAMD,"/mZ=",mZ_TAMD,"/pmf/pmf_T=",T_TAMD,".Zhist_0.2",sep="")

cat(name)

xrange <- c(-0.5,1.0,3)
xrange.axis <- c(-0.5,1.0,3)

yrange <- c(-1.0,0.5,3)
yrange.axis <- c(-1.0,0.5,3)
	
felwFillConMapwrangeworwoaxis(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,yaflag="f",mai=c(0.8,0.0,0.0,0.1),norm="T")

felwFillConBarwmai(name,label.x=label.x,label.y=label.y,level=level,mai=c(7.8,0.5,0.5,0.2))

