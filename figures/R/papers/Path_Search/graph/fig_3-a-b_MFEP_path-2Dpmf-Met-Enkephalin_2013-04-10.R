source("~/calspa/GOLMAA/SH3/R/filled.contour.mod.R")
source("~/calspa/GOLMAA/SH3/R/fel3.R")
source("~/Rspa/plRamachad3.R")

dir <- "~/papers/Path_Search"

title <- NULL

label.x <- expression(paste(d1))
label.y <- expression(paste(d2))

name.path <- paste("/home/yamamori/calspa/MFEP/MetEnk/pmf2MGaussian_2013-04-10_bo_e_CG-FG_CMD_NH_2012-06-04/path_TAA=300_TCG1=370_TCG2=370_1000_0_0_K=12_@0.06,0.15-@0.17,0.06.txt",sep='')
name.out <- paste(dir,"/eps/","fig_3-a-b_MFEP_path-2Dpmf-Met-Enkephalin_2013-04-10",sep='')
   
level <- seq(0.0,5,1)

file.name <- paste(name.out,'.eps',sep='')
postscript(file.name,width=5.6,height=2.4,horizontal=FALSE,onefile=FALSE,paper="special")

xrange <- c(0.0,0.3,3)
yrange <- c(0.0,0.3,3)

par(oma = c(0,0,0,0) )
source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange.R")
source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange_wy.R")

source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMapwrangewxwy_wpath.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMapwrangewxwoy_wpath.R")

source("~/papers/CG-FG_TACCM_REMD/FillConBar.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConBar.R")

nf <- layout(matrix(c(1,2,3),1,3,byrow=TRUE),c(3.7,3,1.5),c(1))

#par(mar =c(5,5,2,2) )
par(cex.axis=1.0)
par(cex.lab=2.0)

name.pmf <- paste("/home/yamamori/calspa/MFEP/MetEnk/pmf2MGaussian_2013-04-10_bo_e_CG-FG_CMD_NH_2012-06-04/pmf2MGaussian_TAA=300_TCG1=370_TCG2=370_1000_0_0_K=12_@0.06,0.15-@0.17,0.06.txt",sep='')
cat(name.pmf," ",name.path,'\n')
felwFillConMapwrangewxwywpath(name.pmf,name.path,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",mar=c(5,5,1,0),xrange=xrange,yrange=yrange,plot.axis="yes")

name.pmf <- paste("/home/yamamori/calspa/MFEP/MetEnk/pmf2MGaussian_2013-04-10_bo_e_CG-FG_CMD_NH_2012-06-04/pmf_TAA=300_TCG1_370_TCG2_370_TZ_1400_0.01_1000_0_0_AA_bo10000ps_2",sep='')
cat(name.pmf," ",name.path,'\n')
felwFillConMapwrangewxwoywpath(name.pmf,name.path,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,mar=c(5,0,1,0),plot.axis="yes")

par(mar=c(5,4,1,1))

felwFillConBar(name.pmf,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

