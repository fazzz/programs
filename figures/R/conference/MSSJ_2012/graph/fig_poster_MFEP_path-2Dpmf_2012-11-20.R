source("~/calspa/GOLMAA/SH3/R/filled.contour.mod.R")
source("~/calspa/GOLMAA/SH3/R/fel3.R")
source("~/Rspa/plRamachad3.R")

dir <- "~/gakkai/MSSJ_2012"

title <- NULL

label.x <- expression(paste(phi))
label.y <- expression(paste(psi))

name.path <- paste("/home/yamamori/calspa/MFEP/AD//pmf2MGaussian_2012-10-10_bo_e_CG-FG_NH_2012-07-27/tau=1.0/mZ=100.00/TZ=700/SB_KZMAX=5000_NR4_woeljd0.001/nEX=1000/fq=10/path_TAA=300_TCG=300_0_5000_CG_K=16_@-1.4,1.1-@1.1,-0.8.txt",sep='')
name.out <- paste(dir,"/tiff/","fig_poster_MFEP_path-2Dpmf_2012-11-20",sep='')
   
level <- seq(0.0,10,1)

file.name <- paste(name.out,'.tiff',sep='')
tiff(file.name,width=700,height=300)

par(oma = c(0,0,0,0) )
source("~/papers/CG-FG_TACCM_REMD/FillConMap.R")
source("~/papers/CG-FG_TACCM_REMD/FillConBar.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap_wpath.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConBar.R")

nf <- layout(matrix(c(1,2,3),1,3,byrow=TRUE),c(3,3,1),c(1))

par(mar =c(5,5,2,2) )
par(cex.axis=2.0)
par(cex.lab=2.0)

name.pmf <- paste("/home/yamamori/calspa/MFEP/AD//pmf2MGaussian_2012-10-10_bo_e_CG-FG_NH_2012-07-27/tau=1.0/mZ=100.00/TZ=700/SB_KZMAX=5000_NR4_woeljd0.001/nEX=1000/fq=10/pmf2MGaussian_TAA=300_TCG=300_0_5000_CG_K=16.txt",sep='')
cat(name.pmf," ",name.path,'\n')
felwFillConMapwpath(name.pmf,name.path,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

name.pmf <- paste("/home/yamamori/calspa/TACCM_CGAAREMD/AD/e_CG-FG_NH_2012-07-27/tau=1.0/mZ=100.00/TZ=700/SB_KZMAX=5000_NR4_woeljd0.001/nEX=1000/fq=10/pmf/pmf_TAA=300_TCG_300_TZ_700_0.3_0_5000_CG_4_2012-08-21",sep='')
cat(name.pmf," ",name.path,'\n')
felwFillConMapwpath(name.pmf,name.path,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

felwFillConBar(name.pmf,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

dev.off()
