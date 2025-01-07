level <- seq(0,5,1)

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

level <- seq(0.0,5,1)

name.out <- paste(dir,"/eps/","fig_2-a-b_2DFEL-MixedGaussians-Met-Enkephalin_2013-04-10",sep='')
file.name <- paste(name.out,'.eps',sep='')
postscript(file.name,width=5.6,height=2.0,horizontal=FALSE,onefile=FALSE,paper="special")

xrange <- c(0.0,0.3,3)
yrange <- c(0.0,0.3,3)

par(oma = c(0,0,0,0) )
source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange.R")
source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange_wy.R")
source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange_wox_wy.R")
source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange_woxyaxis.R")

source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap_wrange.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap_wrange_wy.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap_wrange_woxwy.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap_wrange_woxyaxis.R")

source("~/papers/CG-FG_TACCM_REMD/FillConBar.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConBar.R")

#nf <- layout(matrix(c(1,2,3,4,5,6),1,6,byrow=TRUE),c(3,3,3,3,3,1),c(1))
nf <- layout(matrix(c(1,2,3,4,5),1,5,byrow=TRUE),c(4.5,3,3,3,1.5),c(1))

#par(mar =c(5,5,2,2) )
par(cex.axis=1.0)
par(cex.lab=2.0)
label.x <- expression(paste(d1))
label.y <- expression(paste(d2))

#cat("/home/yamamori/calspa/MFEP/MetEnk/pmf2MGaussian_2013-04-10_bo_e_CG-FG_CMD_NH_2012-06-04/pmf2MGaussian_TAA=300_TCG1=370_TCG2=370_1000_0_0_K=4_@0.06,0.15-@0.17,0.06.txt",'\n')
#felwFillConMap("/home/yamamori/calspa/MFEP/MetEnk/pmf2MGaussian_2013-04-10_bo_e_CG-FG_CMD_NH_2012-06-04/pmf2MGaussian_TAA=300_TCG1=370_TCG2=370_1000_0_0_K=4_@0.06,0.15-@0.17,0.06.txt",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

cat("/home/yamamori/calspa/MFEP/MetEnk/pmf2MGaussian_2013-04-10_bo_e_CG-FG_CMD_NH_2012-06-04/pmf2MGaussian_TAA=300_TCG1=370_TCG2=370_1000_0_0_K=6_@0.06,0.15-@0.17,0.06.txt",'\n')
felwFillConMapwrangewy("/home/yamamori/calspa/MFEP/MetEnk/pmf2MGaussian_2013-04-10_bo_e_CG-FG_CMD_NH_2012-06-04/pmf2MGaussian_TAA=300_TCG1=370_TCG2=370_1000_0_0_K=6_@0.06,0.15-@0.17,0.06.txt",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",mar=c(5,5,2,0),xrange=xrange,yrange=yrange,plot.axis="yes")

cat("/home/yamamori/calspa/MFEP/MetEnk/pmf2MGaussian_2013-04-10_bo_e_CG-FG_CMD_NH_2012-06-04/pmf2MGaussian_TAA=300_TCG1=370_TCG2=370_1000_0_0_K=12_@0.06,0.15-@0.17,0.06.txt",'\n')
felwFillConMapwrange("/home/yamamori/calspa/MFEP/MetEnk/pmf2MGaussian_2013-04-10_bo_e_CG-FG_CMD_NH_2012-06-04/pmf2MGaussian_TAA=300_TCG1=370_TCG2=370_1000_0_0_K=12_@0.06,0.15-@0.17,0.06.txt",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",mar=c(5,0,2,0),xrange=xrange,yrange=yrange,plot.axis="yes")

cat("/home/yamamori/calspa/MFEP/MetEnk/pmf2MGaussian_2013-04-10_bo_e_CG-FG_CMD_NH_2012-06-04/pmf2MGaussian_TAA=300_TCG1=370_TCG2=370_1000_0_0_K=13_@0.06,0.15-@0.17,0.06.txt",'\n')
felwFillConMapwrange("/home/yamamori/calspa/MFEP/MetEnk/pmf2MGaussian_2013-04-10_bo_e_CG-FG_CMD_NH_2012-06-04/pmf2MGaussian_TAA=300_TCG1=370_TCG2=370_1000_0_0_K=13_@0.06,0.15-@0.17,0.06.txt",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",mar=c(5,0,2,0),xrange=xrange,yrange=yrange,plot.axis="yes")

cat("/home/yamamori/calspa/MFEP/MetEnk/pmf2MGaussian_2013-04-10_bo_e_CG-FG_CMD_NH_2012-06-04/pmf2MGaussian_TAA=300_TCG1=370_TCG2=370_1000_0_0_K=14_@0.06,0.15-@0.17,0.06.txt",'\n')

felwFillConMapwrange("/home/yamamori/calspa/MFEP/MetEnk/pmf2MGaussian_2013-04-10_bo_e_CG-FG_CMD_NH_2012-06-04/pmf2MGaussian_TAA=300_TCG1=370_TCG2=370_1000_0_0_K=14_@0.06,0.15-@0.17,0.06.txt",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",mar=c(5,0,2,0),xrange=xrange,yrange=yrange,plot.axis="yes")

par(mar=c(5,2,2,1))

felwFillConBar("/home/yamamori/calspa/MFEP/MetEnk/pmf2MGaussian_2013-04-10_bo_e_CG-FG_CMD_NH_2012-06-04/pmf2MGaussian_TAA=300_TCG1=370_TCG2=370_1000_0_0_K=4_@0.06,0.15-@0.17,0.06.txt",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

