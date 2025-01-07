source("~/calspa/GOLMAA/SH3/R/filled.contour.mod.R")
source("~/calspa/GOLMAA/SH3/R/fel3.R")
source("~/Rspa/plRamachad3.R")

dir <- "~/papers/Path_Search"

title <- NULL

label.x <- expression(paste(phi))
label.y <- expression(paste(psi))

name.path <- paste("/home/yamamori/calspa/MFEP/AD//pmf2MGaussian_2012-10-10_bo_e_CG-FG_NH_2013-03-21/tau=1.0/mZ=50000.00/TZ=750/SB_NR4_wrefd0.1/nEX=1000/fq=10/path_TAA=300_TCG=300_0_2000_CG_K=19_@-1.4,1.1-@1.1,-0.8_wx-2.5-2.5_wy-2.5-2.5.txt",sep='')
name.out <- paste(dir,"/eps/","fig_3-a-b_MFEP_path-2Dpmf_2013-03-21",sep='')
   
level <- seq(0.0,10,1)

file.name <- paste(name.out,'.eps',sep='')
postscript(file.name,width=5.6,height=2.4,horizontal=FALSE,onefile=FALSE,paper="special")

xrange <- c(-2.5,2.5,5)
yrange <- c(-2.5,2.5,5)

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

name.pmf <- paste("/home/yamamori/calspa/MFEP/AD//pmf2MGaussian_2012-10-10_bo_e_CG-FG_NH_2013-03-21/tau=1.0/mZ=50000.00/TZ=750/SB_NR4_wrefd0.1/nEX=1000/fq=10/pmf2MGaussian_TAA=300_TCG=300_0_2000_CG_K=19_wx-2.5-2.5_wy-2.5-2.5.txt",sep='')
cat(name.pmf," ",name.path,'\n')
felwFillConMapwrangewxwywpath(name.pmf,name.path,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",mar=c(5,5,1,0),xrange=xrange,yrange=yrange,plot.axis="yes")

name.pmf <- paste("/home/yamamori/calspa/TACCM_CGAAREMD/AD/e_CG-FG_NH_2012-07-27/tau=1.0/mZ=50000.00/TZ=750/SB_NR4_wrefd0.1/nEX=1000/fq=10/pmf/pmf_TAA=300_TCG_300_TZ_750_0.3_0_2000_CG_4_2012-08-21_wrange_-2.5-2.5_-2.5-2.5",sep='')
cat(name.pmf," ",name.path,'\n')
felwFillConMapwrangewxwoywpath(name.pmf,name.path,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,mar=c(5,0,1,0),plot.axis="yes")

par(mar=c(5,4,1,1))

felwFillConBar(name.pmf,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

