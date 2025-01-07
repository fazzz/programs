fact.x <- 1
fact.y <- 1
fact.p <- 1.0/pi

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

dir <- "~/papers/Path_Search"

title <- NULL

name <- paste("dummy",sep='')

level <- seq(0.0,10,1)

name.out <- paste(dir,"/eps/","fig_2-a-b-c-d_2DFEL-MixedGaussians_2013-06-14",sep='')
file.name <- paste(name.out,'.eps',sep='')
postscript(file.name,width=6.4,height=2.0,horizontal=FALSE,onefile=FALSE,paper="special")

xrange <- c(-1.0,1.0,5)
xrange.axis <- c(-1.0,1.0,5)
yrange <- c(-1.0,1.0,5)
yrange.axis <- c(-1.0,1.0,5)

par(oma = c(0,0,0,0) )
source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange.R")
source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange_wy.R")
source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange_wox_wy.R")
source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange_woxyaxis.R")

source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap_wrange.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap_wrange_wy.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap_wrange_woxwy.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap_wrange_woxyaxis.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMapwrangewxwy_wpath.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMapwrangewxwoy_wpath.R")

source("~/papers/CG-FG_TACCM_REMD/FillConBar.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConBar.R")

nf <- layout(matrix(c(1,2,3,4,5),1,5,byrow=TRUE),c(4.5,3,3,3,1.5),c(1))

#par(mar =c(5,5,2,2) )
par(cex.axis=1.44)
par(cex.lab=1.44)
#label.x <- expression(paste(phi))
#label.y <- expression(paste(psi))

label.x <- " " # expression(paste(phi,"(radian/",pi,")"))
label.y <- expression(paste(psi,"(radian/",pi,")"))

name.pmf <- paste("/home/yamamori/calspa/TACCM_CGAAREMD/AD/e_CG-FG_NH_2012-07-27/tau=1.0/mZ=100.00/TZ=750/SB_KZMAX=5000_NR8-2_woeljd0.001/nEX=1000/fq=10/pmf/pmf_TAA=300_TCG_300_TZ_750_0.3_0_5000_CG_4_2012-08-21",sep='')
name.path <- paste("/home/yamamori/calspa/MFEP/AD/pmf2MGaussian_2012-10-10_bo_e_CG-FG_NH_2013-05-04/tau=1.0/mZ=100.00/TZ=750/SB_KZMAX=5000_NR8-2_woeljd0.001/nEX=1000/fq=10/path_TAA=300_TCG=300_0_5000_CG_K=18_@-1.4,1.1-@1.1,-0.8_wx-4.0-4.0_wy-4.0-4.0.txt",sep='')

xrange <- c(-1.0,1.0,4)
xrange.axis <- c(-1.0,1.0,4)
yrange <- c(-1.0,1.0,4)
yrange.axis <- c(-1.0,1.0,4)

cat("/home/yamamori/calspa/TACCM_CGAAREMD/AD/e_CG-FG_NH_2012-07-27/tau=1.0/mZ=100.00/TZ=750/SB_KZMAX=5000_NR8-2_woeljd0.001/nEX=1000/fq=10/pmf/pmf_TAA=300_TCG_300_TZ_750_0.3_0_5000_CG_4_2012-08-21",'\n')
felwFillConMapwrangewxwywpath(name.pmf,name.path,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",mar=c(5,5,1,0),xrange=xrange,yrange=yrange,plot.axis="yes",norm="T")

xrange <- c(-0.5,1.0,3)
xrange.axis <- c(-0.5,1.0,3)
yrange <- c(-1.0,1.0,4)
yrange.axis <- c(-1.0,1.0,5)

cat("/home/yamamori/calspa/MFEP/AD/pmf2MGaussian_2012-10-10_bo_e_CG-FG_NH_2013-05-04/tau=1.0/mZ=100.00/TZ=750/SB_KZMAX=5000_NR8-2_woeljd0.001/nEX=1000/fq=10/pmf2MGaussian_TAA=300_TCG=300_0_5000_CG_K=6_wx-4.0-4.0_wy-4.0-4.0_dx=0.01_dy=0.01_4.txt",'\n')
felwFillConMapwrange("/home/yamamori/calspa/MFEP/AD/pmf2MGaussian_2012-10-10_bo_e_CG-FG_NH_2013-05-04/tau=1.0/mZ=100.00/TZ=750/SB_KZMAX=5000_NR8-2_woeljd0.001/nEX=1000/fq=10/pmf2MGaussian_TAA=300_TCG=300_0_5000_CG_K=6_wx-4.0-4.0_wy-4.0-4.0_dx=0.01_dy=0.01_4.txt",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",mar=c(5,0,1,0),xrange=xrange,yrange=yrange,plot.axis="yes",norm="T")

xrange <- c(-0.5,1.0,3)
xrange.axis <- c(-0.5,1.0,3)
yrange <- c(-1.0,1.0,4)
yrange.axis <- c(-1.0,1.0,5)

cat("/home/yamamori/calspa/MFEP/AD/pmf2MGaussian_2012-10-10_bo_e_CG-FG_NH_2013-05-04/tau=1.0/mZ=100.00/TZ=750/SB_KZMAX=5000_NR8-2_woeljd0.001/nEX=1000/fq=10/pmf2MGaussian_TAA=300_TCG=300_0_5000_CG_K=17_wx-4.0-4.0_wy-4.0-4.0_dx=0.01_dy=0.01_4.txt",'\n')
felwFillConMapwrange("/home/yamamori/calspa/MFEP/AD/pmf2MGaussian_2012-10-10_bo_e_CG-FG_NH_2013-05-04/tau=1.0/mZ=100.00/TZ=750/SB_KZMAX=5000_NR8-2_woeljd0.001/nEX=1000/fq=10/pmf2MGaussian_TAA=300_TCG=300_0_5000_CG_K=17_wx-4.0-4.0_wy-4.0-4.0_dx=0.01_dy=0.01_4.txt",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",mar=c(5,0,1,0),xrange=xrange,yrange=yrange,plot.axis="yes",norm="T")

xrange <- c(-0.5,1.0,3)
xrange.axis <- c(-0.5,1.0,3)
yrange <- c(-1.0,1.0,4)
yrange.axis <- c(-1.0,1.0,5)

cat("/home/yamamori/calspa/MFEP/AD/pmf2MGaussian_2012-10-10_bo_e_CG-FG_NH_2013-05-04/tau=1.0/mZ=100.00/TZ=750/SB_KZMAX=5000_NR8-2_woeljd0.001/nEX=1000/fq=10/pmf2MGaussian_TAA=300_TCG=300_0_5000_CG_K=18_wx-4.0-4.0_wy-4.0-4.0_dx=0.01_dy=0.01_4.txt",'\n')
felwFillConMapwrange("/home/yamamori/calspa/MFEP/AD/pmf2MGaussian_2012-10-10_bo_e_CG-FG_NH_2013-05-04/tau=1.0/mZ=100.00/TZ=750/SB_KZMAX=5000_NR8-2_woeljd0.001/nEX=1000/fq=10/pmf2MGaussian_TAA=300_TCG=300_0_5000_CG_K=18_wx-4.0-4.0_wy-4.0-4.0_dx=0.01_dy=0.01_4.txt",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",mar=c(5,0,1,0),xrange=xrange,yrange=yrange,plot.axis="yes",norm="T")

par(mar=c(5,2,1,1))

felwFillConBar("/home/yamamori/calspa/TACCM_CGAAREMD/AD/e_CG-FG_NH_2012-07-27/tau=1.0/mZ=100.00/TZ=750/SB_KZMAX=5000_NR8-2_woeljd0.001/nEX=1000/fq=10/pmf/pmf_TAA=300_TCG_300_TZ_750_0.3_0_5000_CG_4_2012-08-21",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

