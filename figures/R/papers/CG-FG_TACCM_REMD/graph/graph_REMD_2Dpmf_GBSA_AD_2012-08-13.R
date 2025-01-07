
source("~/calspa/GOLMAA/SH3/R/filled.contour.mod.R")
source("~/calspa/GOLMAA/SH3/R/fel3.R")
source("~/Rspa/plRamachad3.R")

dir0 <- "~/calspa/refcalc/REMD/AD"
dirbase <- paste(dir0,"/s_REGB_2012-07-23_",ff,sep="")

nTZs<-1

title<-NULL

label.x <- expression(paste(phi))
label.y <- expression(paste(psi))

name <- paste(dirbase,"/",pname,"/nEX=",numEX,"/freq=",TLbase,"ps","/pmf/pmf_ADg",sep="")

cat(name,'\n')

file.name <- paste(name.out,'.tiff',sep='')
tiff(file.name,width=390,height=320)
fel3(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")
dev.off()
