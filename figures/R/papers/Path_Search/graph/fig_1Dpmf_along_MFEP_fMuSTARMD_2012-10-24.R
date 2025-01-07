source("~/calspa/GOLMAA/SH3/R/filled.contour.mod.R")
source("~/calspa/GOLMAA/SH3/R/fel3.R")
source("~/Rspa/plRamachad3.R")

dir <- "~/papers/Path_Search"

title <- NULL

label.x <- expression(paste(phi))
label.y <- expression(paste(psi))

name.path <- paste("/home/yamamori/calspa/MFEP/AD//pmf2MGaussian_2012-10-10_bo_e_CG-FG_NH_2012-07-27/tau=1.0/mZ=50000.00/TZ=750/SB_NR4_wrefd0.1/nEX=1000/fq=10/path_TAA=300_TCG=300_0_2000_CG_K=17_@-1.4,1.1-@1.1,-0.8.txt",sep='')
name.out <- paste(dir,"/tiff/","fig_1Dpmf_along_MFEP_fMuSTARMD_2012-10-24",sep='')
   
level <- seq(0.0,10,1)

file.name <- paste(name.out,'.tiff',sep='')
tiff(file.name,width=700,height=300)

par(oma = c(0,0,0,0) )
source("~/papers/CG-FG_TACCM_REMD/FillConMap.R")
source("~/papers/CG-FG_TACCM_REMD/FillConBar.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap_wpath.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConBar.R")

nf <- layout(matrix(c(1,1,2,3),2,2,byrow=TRUE),c(3,1),c(1,1))

par(mar =c(5,5,2,2) )
par(cex.axis=2.0)
par(cex.lab=2.0)

label.x=" "
label.y="pmf"

xrange <- c(0,1.0,10)
xrange.axis <- xrange
yrange <- c(0,10,10)
yrange.axis <- yrange

fact.x <- 1
fact.y <- 1
ave.x <- 1

id.xs <- c(1,1)
id.ys <- c(2,2)
ids.ys <- c(3,3)
iro <- c(1,2)
senshu <- c(1,1)
tenshu <- c(3,3)
hutosa <- c(1,1)
is.leg <- 0

source("~/Rspa/plmGeneral_wsrange.R")

########################################################
# plmGeneralwsrange(data.names="/home/yamamori/calspa/MFEP/AD//pmf2MGaussian_2012-10-10_bo_e_CG-FG_NH_2012-07-27/tau=1.0/mZ=50000.00/TZ=750/SB_NR4_wrefd0.1/nEX=1000/fq=10/path2_TAA=300_TCG=300_0_2000_CG_K=17_@-1.4,1.1-@1.1,-0.8.txt",     #
#                   sd.names="/home/yamamori/calspa/MFEP/AD//pmf2MGaussian_2012-10-10_bo_e_CG-FG_NH_2012-07-27/tau=1.0/mZ=50000.00/TZ=750/SB_NR4_wrefd0.1/nEX=1000/fq=10/path2_TAA=300_TCG=300_0_2000_CG_K=17_@-1.4,1.1-@1.1,-0.8.txt",       #
#                   id.ys=id.ys,		       #
#                   ids.ys=ids.ys,		       #
#                   is.sen=rep(0,20),		       #
#                   label.size=0.5,axis.size=2.0,      #
#                   iro=iro,axis.ft="F",is.header="T", #
#                   sdiro=iro,			       #
#                   xrange=xrange,yrange=yrange,       #
#                   sdyrange=yrange,		       #
#                   warrow="T")			       #
########################################################
#box(lwd=2.0)
  
#axis(2,yaxp=yrange.axis,lwd=2.0,cex.axis=1.5)
#mtext(outer=T,label.y,side=2,line=4.0,cex=1.5)
  
#axis(1,xaxp=xrange.axis,lwd=2.0,cex.axis=1.5)
#mtext(outer=T,label.x,side=1,line=4.0,cex=1.5)

name.pmf <- paste("/home/yamamori/calspa/TACCM_CGAAREMD/AD/e_CG-FG_NH_2012-07-27/tau=1.0/mZ=50000.00/TZ=750/SB_NR4_wrefd0.1/nEX=1000/fq=10/pmf/pmf_TAA=300_TCG_300_TZ_750_0.3_0_2000_CG_4_2012-08-21",sep='')
cat(name.pmf," ",name.path,'\n')
felwFillConMapwpath(name.pmf,name.path,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

name.pmf <- paste("/home/yamamori/calspa/TACCM_CGAAREMD/AD/e_CG-FG_NH_2012-07-27/tau=1.0/mZ=50000.00/TZ=750/SB_NR4_wrefd0.1/nEX=1000/fq=10/pmf/pmf_TAA=300_TCG_300_TZ_750_0.3_0_2000_CG_4_2012-08-21",sep='')
cat(name.pmf," ",name.path,'\n')
felwFillConMapwpath(name.pmf,name.path,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

felwFillConBar(name.pmf,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")

dev.off()
