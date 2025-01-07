fact.x <-180/pi
fact.y <-180/pi
fact.p <-180/pi

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

dir <- "/home/yamamori/defense"

title <- NULL

name <- paste("dummy",sep='')

level <- seq(0.0,10,1)

name.out <- paste(dir,"/eps/","fig_PathSearch_2DFEL-MixedGaussians_wpath_2013-07-01",sep='')
file.name <- paste(name.out,'.eps',sep='')
#postscript(file.name,width=5.8024,height=1.6236,horizontal=FALSE,onefile=FALSE,paper="special")
postscript(file.name,width=6.96228,height=1.94832,horizontal=FALSE,onefile=FALSE,paper="special")

xrange <- c(-180,180,6)
xrange.axis <- c(-180,180,6)
yrange <- c(-180,180,6)
yrange.axis <- c(-180,180,6)

par(oma = c(0,0,0,0) )
source("~/defense/fig/FillConMap.R")
source("~/defense/fig/FillConBar.R")
source("~/defense/fig/fel_FillConMap.R")
source("~/defense/fig/fel_FillConMap_wpath.R")
source("~/defense/fig/fel_FillConBar.R")

nf <- layout(matrix(c(1,2,3,4,5),1,5,byrow=TRUE),c(1.5,1.1,1.1,1.1,1.0),c(1))

#par(mar =c(5,5,2,2) )
par(cex.axis=1.0)
par(cex.lab=2.0)
label.x <- ""
label.y <- ""

cat("/home/yamamori/calspa/TACCM_CGAAREMD/AD/e_CG-FG_NH_2012-07-27/tau=1.0/mZ=50000.00/TZ=750/SB_NR4_wrefd0.1/nEX=1000/fq=10/pmf/pmf_TAA=300_TCG_300_TZ_750_0.3_0_2000_CG_4_2012-08-21_wrange_-2.5-2.5_-2.5-2.5",'\n')
cat("/home/yamamori/calspa/MFEP/AD/pmf2MGaussian_2012-10-10_bo_e_CG-FG_NH_2013-03-21/tau=1.0/mZ=50000.00/TZ=750/SB_NR4_wrefd0.1/nEX=1000/fq=10/path_TAA=300_TCG=300_0_2000_CG_K=19_@-1.4,1.1-@1.1,-0.8_wx-2.5-2.5_wy-2.5-2.5.txt",'\n')

felwFillConMapwpath("/home/yamamori/calspa/TACCM_CGAAREMD/AD/e_CG-FG_NH_2012-07-27/tau=1.0/mZ=50000.00/TZ=750/SB_NR4_wrefd0.1/nEX=1000/fq=10/pmf/pmf_TAA=300_TCG_300_TZ_750_0.3_0_2000_CG_4_2012-08-21_wrange_-2.5-2.5_-2.5-2.5","/home/yamamori/calspa/MFEP/AD/pmf2MGaussian_2012-10-10_bo_e_CG-FG_NH_2013-03-21/tau=1.0/mZ=50000.00/TZ=750/SB_NR4_wrefd0.1/nEX=1000/fq=10/path_TAA=300_TCG=300_0_2000_CG_K=19_@-1.4,1.1-@1.1,-0.8_wx-2.5-2.5_wy-2.5-2.5.txt",label.x=label.x,label.y=label.y,level=level,
               title=title,kt="F/kBT",
               xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,
               mai=c(0.5,0.5,0.1,0))

cat("/home/yamamori/calspa/MFEP/AD/pmf2MGaussian_2012-10-10_bo_e_CG-FG_NH_2013-03-21/tau=1.0/mZ=50000.00/TZ=750/SB_NR4_wrefd0.1/nEX=1000/fq=10/pmf2MGaussian_TAA=300_TCG=300_0_2000_CG_K=4_wx-2.5-2.5_wy-2.5-2.5.txt",'\n')
felwFillConMap("/home/yamamori/calspa/MFEP/AD/pmf2MGaussian_2012-10-10_bo_e_CG-FG_NH_2013-03-21/tau=1.0/mZ=50000.00/TZ=750/SB_NR4_wrefd0.1/nEX=1000/fq=10/pmf2MGaussian_TAA=300_TCG=300_0_2000_CG_K=4_wx-2.5-2.5_wy-2.5-2.5.txt",label.x=label.x,label.y=label.y,level=level,
               title=title,kt="F/kBT",
               xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,
               yaflag="F",
               mai=c(0.5,0.0,0.1,0))

cat("/home/yamamori/calspa/MFEP/AD/pmf2MGaussian_2012-10-10_bo_e_CG-FG_NH_2013-03-21/tau=1.0/mZ=50000.00/TZ=750/SB_NR4_wrefd0.1/nEX=1000/fq=10/pmf2MGaussian_TAA=300_TCG=300_0_2000_CG_K=10_wx-2.5-2.5_wy-2.5-2.5.txt",'\n')
felwFillConMap("/home/yamamori/calspa/MFEP/AD/pmf2MGaussian_2012-10-10_bo_e_CG-FG_NH_2013-03-21/tau=1.0/mZ=50000.00/TZ=750/SB_NR4_wrefd0.1/nEX=1000/fq=10/pmf2MGaussian_TAA=300_TCG=300_0_2000_CG_K=10_wx-2.5-2.5_wy-2.5-2.5.txt",label.x=label.x,label.y=label.y,level=level,
               title=title,kt="F/kBT",
               xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,
               yaflag="F",
               mai=c(0.5,0.0,0.1,0))

cat("/home/yamamori/calspa/MFEP/AD/pmf2MGaussian_2012-10-10_bo_e_CG-FG_NH_2013-03-21/tau=1.0/mZ=50000.00/TZ=750/SB_NR4_wrefd0.1/nEX=1000/fq=10/pmf2MGaussian_TAA=300_TCG=300_0_2000_CG_K=18_wx-2.5-2.5_wy-2.5-2.5.txt",'\n')
felwFillConMap("/home/yamamori/calspa/MFEP/AD/pmf2MGaussian_2012-10-10_bo_e_CG-FG_NH_2013-03-21/tau=1.0/mZ=50000.00/TZ=750/SB_NR4_wrefd0.1/nEX=1000/fq=10/pmf2MGaussian_TAA=300_TCG=300_0_2000_CG_K=18_wx-2.5-2.5_wy-2.5-2.5.txt",label.x=label.x,label.y=label.y,level=level,
               title=title,kt="F/kBT",
               xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,
               yaflag="F",
               mai=c(0.5,0.0,0.1,0))

felwFillConBar("/home/yamamori/calspa/TACCM_CGAAREMD/AD/e_CG-FG_NH_2012-07-27/tau=1.0/mZ=50000.00/TZ=750/SB_NR4_wrefd0.1/nEX=1000/fq=10/pmf/pmf_TAA=300_TCG_300_TZ_750_0.3_0_2000_CG_4_2012-08-21_wrange_-2.5-2.5_-2.5-2.5",label.x=label.x,label.y=label.y,level=level,
               mai=c(0.94,0.75,0.1,0.5))

