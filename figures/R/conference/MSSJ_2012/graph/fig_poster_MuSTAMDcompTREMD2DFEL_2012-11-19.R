source("~/calspa/GOLMAA/SH3/R/filled.contour.mod.R")
source("~/calspa/GOLMAA/SH3/R/fel3.R")
source("~/Rspa/plRamachad3.R")

dir <- "~/gakkai/MSSJ_2012"

title <- NULL

label.x <- expression(paste(phi))
label.y <- expression(paste(psi))

name.path <- paste("",sep='')
name.out <- paste(dir,"/tiff/","fig_poster_MuSTAMDcompTREMD2DFEL_2012-11-19",sep='')
   
level <- seq(0.0,20,2.4)

file.name <- paste(name.out,'.tiff',sep='')
tiff(file.name,width=720,height=300)

par(oma = c(0,0,0,0) )
source("~/papers/CG-FG_TACCM_REMD/FillConMap.R")
source("~/papers/CG-FG_TACCM_REMD/FillConBar.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConBar.R")

nf <- layout(matrix(c(1,2,3),1,3,byrow=TRUE),c(3,3,1),c(1))

#par(mar =c(5,5,2,2) )
#par(cex.axis=2.0)
#par(cex.lab=2.0)

par(mar =c(4,4,2,2) )
par(cex.axis=2.0)
par(cex.lab=2.0)

felwFillConMap("/home/yamamori/calspa/refcalc/REMD/AD/s_REVAC_2012-07-23_ff99SB/f300t400/nEX=100000/freq=1ps/pmf/pmf_ADv",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

name.pmf <- paste("/home/yamamori/calspa/TACCM_CGAAREMD/AD/e_CG-FG_NH_2012-07-27/tau=1.0/mZ=100.00/TZ=700/SB_KZMAX=5000_NR4_woeljd0.001/nEX=1000/fq=10/pmf/pmf_TAA=300_TCG_300_TZ_700_0.3_0_5000_CG_4_2012-08-21",sep='')
felwFillConMap(name.pmf,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

felwFillConBar(name.pmf,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

dev.off()
