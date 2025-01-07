
dir <- "~/defense/"

title <- NULL

label.x <- expression(paste(phi))
label.y <- expression(paste(psi))

name.out <- paste(dir,"/eps/","fig_MFEP_pathon2DFEL_2013-07-01",sep='')
   
level <- seq(0.0,10,1)

file.name <- paste(name.out,'.eps',sep='')
postscript(file.name,width=4.1496,height=2.6472,horizontal=FALSE,onefile=FALSE,paper="special")

xrange <- c(-180,180,6)
xrange.axis <- c(-180,180,6)
yrange <- c(-180,180,6)
yrange.axis <- c(-180,180,6)

par(oma = c(0,0,0,0) )
source("~/defense/fig/FillConMap.R")
source("~/defense/fig/FillConBar.R")

source("~/defense/fig/fel_FillConMap_wpath.R")
source("~/defense/fig/fel_FillConBar.R")

nf <- layout(matrix(c(1,2),1,2,byrow=TRUE),c(6.47,3.07),c(1))

par(cex.axis=1.0)
par(cex.lab=1.5)

fact.x <- 180/pi
fact.y <- 180/pi
fact.p <- 180/pi

felwFillConMapwpath("/home/yamamori/calspa/TACCM_CGAAREMD/AD/e_CG-FG_NH_2012-07-27/tau=1.0/mZ=50000.00/TZ=750/SB_NR4_wrefd0.1/nEX=1000/fq=10/pmf/pmf_TAA=300_TCG_300_TZ_750_0.3_0_2000_CG_4_2012-08-21","/home/yamamori/calspa/MFEP/AD//pmf2MGaussian_2012-10-10_bo_e_CG-FG_NH_2013-03-21/tau=1.0/mZ=50000.00/TZ=750/SB_NR4_wrefd0.1/nEX=1000/fq=10/path_TAA=300_TCG=300_0_2000_CG_K=_@-1.4,1.1-@1.1,-0.8_wx-4.0-4.0_wy-4.0-4.0.txt",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",mai=c(0.5,0.5,0.1,0.0),xrange=xrange,yrange=yrange,plot.axis="yes")

felwFillConBar("/home/yamamori/calspa/TACCM_CGAAREMD/AD/e_CG-FG_NH_2012-07-27/tau=1.0/mZ=50000.00/TZ=750/SB_NR4_wrefd0.1/nEX=1000/fq=10/pmf/pmf_TAA=300_TCG_300_TZ_750_0.3_0_2000_CG_4_2012-08-21",label.x=label.x,label.y=label.y,level=level,mai=c(1.348,0.75,0.1,0.75))


dir <- "~/defense/"

title <- NULL

label.x <- expression(paste(phi))
label.y <- expression(paste(psi))

name.out <- paste(dir,"/eps/","fig_MFEP_pathon2DFEL_2013-07-01",sep='')
   
level <- seq(0.0,10,1)

file.name <- paste(name.out,'.eps',sep='')
postscript(file.name,width=4.1496,height=2.6472,horizontal=FALSE,onefile=FALSE,paper="special")

xrange <- c(-180,180,6)
xrange.axis <- c(-180,180,6)
yrange <- c(-180,180,6)
yrange.axis <- c(-180,180,6)

par(oma = c(0,0,0,0) )
source("~/defense/fig/FillConMap.R")
source("~/defense/fig/FillConBar.R")

source("~/defense/fig/fel_FillConMap_wpath.R")
source("~/defense/fig/fel_FillConBar.R")

nf <- layout(matrix(c(1,2),1,2,byrow=TRUE),c(6.47,3.07),c(1))

par(cex.axis=1.0)
par(cex.lab=1.5)

fact.x <- 180/pi
fact.y <- 180/pi
fact.p <- 180/pi

felwFillConMapwpath("/home/yamamori/calspa/TACCM_CGAAREMD/AD/e_CG-FG_NH_2012-07-27/tau=1.0/mZ=50000.00/TZ=750/SB_NR4_wrefd0.1/nEX=1000/fq=10/pmf/pmf_TAA=300_TCG_300_TZ_750_0.3_0_2000_CG_4_2012-08-21","/home/yamamori/calspa/MFEP/AD//pmf2MGaussian_2012-10-10_bo_e_CG-FG_NH_2013-03-21/tau=1.0/mZ=50000.00/TZ=750/SB_NR4_wrefd0.1/nEX=1000/fq=10/path_TAA=300_TCG=300_0_2000_CG_K=_@-1.4,1.1-@1.1,-0.8_wx-4.0-4.0_wy-4.0-4.0.txt",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",mai=c(0.5,0.5,0.1,0.0),xrange=xrange,yrange=yrange,plot.axis="yes")

felwFillConBar("/home/yamamori/calspa/TACCM_CGAAREMD/AD/e_CG-FG_NH_2012-07-27/tau=1.0/mZ=50000.00/TZ=750/SB_NR4_wrefd0.1/nEX=1000/fq=10/pmf/pmf_TAA=300_TCG_300_TZ_750_0.3_0_2000_CG_4_2012-08-21",label.x=label.x,label.y=label.y,level=level,mai=c(1.348,0.75,0.1,0.75))


dir <- "~/defense/"

title <- NULL

label.x <- expression(paste(phi))
label.y <- expression(paste(psi))

name.out <- paste(dir,"/eps/","fig_MFEP_pathon2DFEL_2013-07-01",sep='')
   
level <- seq(0.0,10,1)

file.name <- paste(name.out,'.eps',sep='')
postscript(file.name,width=4.1496,height=2.6472,horizontal=FALSE,onefile=FALSE,paper="special")

xrange <- c(-180,180,6)
xrange.axis <- c(-180,180,6)
yrange <- c(-180,180,6)
yrange.axis <- c(-180,180,6)

par(oma = c(0,0,0,0) )
source("~/defense/fig/FillConMap.R")
source("~/defense/fig/FillConBar.R")

source("~/defense/fig/fel_FillConMap_wpath.R")
source("~/defense/fig/fel_FillConBar.R")

nf <- layout(matrix(c(1,2),1,2,byrow=TRUE),c(6.47,3.07),c(1))

par(cex.axis=1.0)
par(cex.lab=1.5)

fact.x <- 180/pi
fact.y <- 180/pi
fact.p <- 180/pi

felwFillConMapwpath("/home/yamamori/calspa/TACCM_CGAAREMD/AD/e_CG-FG_NH_2012-07-27/tau=1.0/mZ=50000.00/TZ=750/SB_NR4_wrefd0.1/nEX=1000/fq=10/pmf/pmf_TAA=300_TCG_300_TZ_750_0.3_0_2000_CG_4_2012-08-21","/home/yamamori/calspa/MFEP/AD//pmf2MGaussian_2012-10-10_bo_e_CG-FG_NH_2013-03-21/tau=1.0/mZ=50000.00/TZ=750/SB_NR4_wrefd0.1/nEX=1000/fq=10/path_TAA=300_TCG=300_0_2000_CG_K=_@-1.4,1.1-@1.1,-0.8_wx-4.0-4.0_wy-4.0-4.0.txt",label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",mai=c(0.5,0.5,0.1,0.0),xrange=xrange,yrange=yrange,plot.axis="yes")

felwFillConBar("/home/yamamori/calspa/TACCM_CGAAREMD/AD/e_CG-FG_NH_2012-07-27/tau=1.0/mZ=50000.00/TZ=750/SB_NR4_wrefd0.1/nEX=1000/fq=10/pmf/pmf_TAA=300_TCG_300_TZ_750_0.3_0_2000_CG_4_2012-08-21",label.x=label.x,label.y=label.y,level=level,mai=c(1.348,0.75,0.1,0.75))

