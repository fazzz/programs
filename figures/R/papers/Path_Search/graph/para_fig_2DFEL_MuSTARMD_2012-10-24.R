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

name.out <- paste(dir,"/tiff/","fig_2DFEL_MuSTARMD_2012-10-24",sep='')
file.name <- paste(name.out,'.tiff',sep='')
tiff(file.name,width=400,height=300)

par(oma = c(0,0,0,0) )
source("~/papers/CG-FG_TACCM_REMD/FillConMap.R")
source("~/papers/CG-FG_TACCM_REMD/FillConBar.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConBar.R")
						
nf <- layout(matrix(c(1,2),1,2,byrow=TRUE),c(3,1),c(1))

par(mar =c(5,5,2,2) )
par(cex.axis=2.0)
par(cex.lab=2.0)
label.x <- expression(paste(phi))
label.y <- expression(paste(psi))

cat("/home/yamamori/calspa/TACCM_CGAAREMD/AD/e_CG-FG_NH_2012-07-27/tau=1.0/mZ=50000.00/TZ=750/SB_NR4_wrefd0.1/nEX=1000/fq=10/pmf/pmf_TAA=300_TCG_300_TZ_750_0.3_0_2000_CG_4_2012-08-21",'\n')
felwFillConMap("/home/yamamori/calspa/TACCM_CGAAREMD/AD/e_CG-FG_NH_2012-07-27/tau=1.0/mZ=50000.00/TZ=750/SB_NR4_wrefd0.1/nEX=1000/fq=10/pmf/pmf_TAA=300_TCG_300_TZ_750_0.3_0_2000_CG_4_2012-08-21",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

felwFillConBar("/home/yamamori/calspa/TACCM_CGAAREMD/AD/e_CG-FG_NH_2012-07-27/tau=1.0/mZ=50000.00/TZ=750/SB_NR4_wrefd0.1/nEX=1000/fq=10/pmf/pmf_TAA=300_TCG_300_TZ_750_0.3_0_2000_CG_4_2012-08-21",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

dev.off()
