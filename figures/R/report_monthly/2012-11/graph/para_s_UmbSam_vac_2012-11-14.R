dir <- "/home/yamamori/calspa/refcalc/UmbSam/AD"
dirout <- "/home/yamamori/Report/2012-11"

title <- NULL

label.x <- expression(paste(phi))
label.y <- expression(paste(psi))

name.out <- paste(dirout,"/tiff/","fig_s_UmbSam_vac_2012-11-14",sep='')
   
level <- seq(0.0,20,2)

file.name <- paste(name.out,'.tiff',sep='')
tiff(file.name,width=420,height=300)

par(oma = c(0,0,0,0) )
source("~/papers/CG-FG_TACCM_REMD/FillConMap.R")
source("~/papers/CG-FG_TACCM_REMD/FillConBar.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConBar.R")

nf <- layout(matrix(c(1,2),1,2,byrow=TRUE),c(3,1),c(1))

par(mar =c(5,5,2,2) )
par(cex.axis=2.0)
par(cex.lab=2.0)

name.pmf <- paste("/home/yamamori/calspa/refcalc/UmbSam/AD/s_UmbSam_vac_2012-11-12_ff99SB/Umb_Nbin=12x12_K=10/pmf/pmf_UmbMD_vac_10ns.txt",sep='')
felwFillConMap(name.pmf,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")
felwFillConBar(name.pmf,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

dev.off()
