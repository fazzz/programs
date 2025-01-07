source("~/calspa/GOLMAA/SH3/R/filled.contour.mod.R")
source("~/calspa/GOLMAA/SH3/R/fel3.R")
source("~/Rspa/plRamachad3.R")

dir <- "~/gakkai/MSSJ_2012"

title <- NULL

label.x <- expression(paste(phi))
label.y <- expression(paste(psi))

name.path <- paste("/home/yamamori/calspa/MFEP/AD//pmf2MGaussian_2012-10-10_bo_e_CG-FG_NH_2012-07-27/tau=1.0/mZ=50000.00/TZ=750/SB_NR4_wrefd0.1/nEX=1000/fq=10/path_TAA=300_TCG=300_0_2000_CG_K=17_@-1.4,1.1-@1.1,-0.8.txt",sep='')
name.out <- paste(dir,"/tiff/","fig_abst_MFEPon2DFEL_2012-10-22",sep='')
   
level <- seq(0.0,10,1)

file.name <- paste(name.out,'.tiff',sep='')
tiff(file.name,width=400,height=300)

par(oma = c(0,0,0,0) )
source("~/papers/CG-FG_TACCM_REMD/FillConMap.R")
source("~/papers/CG-FG_TACCM_REMD/FillConBar.R")
source("~/gakkai/MSSJ_2012/graph/fel_wFillConMap_wpath.R")
source("~/gakkai/MSSJ_2012/graph/fel_wFillConBar.R")

nf <- layout(matrix(c(1,2),1,2,byrow=TRUE),c(3,1),c(1))

par(mar =c(5,5,2,2) )
par(cex.axis=2.0)
par(cex.lab=2.0)

name.pmf <- paste("/home/yamamori/calspa/TACCM_CGAAREMD/AD/e_CG-FG_NH_2012-07-27/tau=1.0/mZ=50000.00/TZ=750/SB_NR4_wrefd0.1/nEX=1000/fq=10/pmf/pmf_TAA=300_TCG_300_TZ_750_0.3_0_2000_CG_4_2012-08-21",sep='')
cat(name.pmf," ",name.path,'\n')
felwFillConMapwpath(name.pmf,name.path,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

felwFillConBar(name.pmf,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

dev.off()
