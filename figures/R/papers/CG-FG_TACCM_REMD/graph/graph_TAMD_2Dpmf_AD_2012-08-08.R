
source("~/calspa/GOLMAA/SH3/R/filled.contour.mod.R")
source("~/calspa/GOLMAA/SH3/R/fel3.R")
source("~/Rspa/plRamachad3.R")

dir0 <- "/home/yamamori/calspa/TACCM/AD"
dirbase <- paste(dir0,"/e_TACCM_NH_2012-07-31_",ff,sep="")

title<-NULL

label.x <- expression(paste(phi))
label.y <- expression(paste(psi))

name <- paste(dirbase,"/tau=",tau,"/TB=",TB,"/KZ=",KZ,"/mZ=",mZ,"/pmf/pmf_T=",T,".Zhist",sep="")

cat(name,'\n')

file.name <- paste(name.out,'.tiff',sep='')
tiff(file.name,width=390,height=320)
#tiff(file.name,width=312,height=256)
#tiff(file.name,width=293,height=240)
fel3(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")
dev.off()

