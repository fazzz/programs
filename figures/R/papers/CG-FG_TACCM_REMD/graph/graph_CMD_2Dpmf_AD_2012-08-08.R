 
source("~/calspa/GOLMAA/SH3/R/filled.contour.mod.R")
source("~/calspa/GOLMAA/SH3/R/fel3.R")
source("~/Rspa/plRamachad3.R")

dir0 <- "~/calspa/refcalc/CMD/AD"
dirbase <- paste(dir0,"/s_CVAC_2012-08-08_",ff,sep="")

nTZs<-1

label.x <- expression(paste(phi))
label.y <- expression(paste(psi))

name <- paste(dirbase,"/anl/pmf_ADv_T",T,"_",width,sep="")

cat(name,'\n')

file.name <- paste(name.out,'.tiff',sep='')
tiff(file.name,width=390,height=320)
fel3(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")
dev.off()
