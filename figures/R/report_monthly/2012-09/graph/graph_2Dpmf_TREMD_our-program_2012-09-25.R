
source("~/calspa/GOLMAA/SH3/R/filled.contour.mod.R")
source("~/calspa/GOLMAA/SH3/R/fel3.R")
source("~/Rspa/plRamachad3.R")

dir0 <- "~/calspa/refcalc/REMD/AD"
dirbase <- paste(dir0,"/s_REVAC_2012-07-23_",ff,sep="")

nTZs<-1

name.title <- ""

title=name.title

label.x=expression(phi)
label.y=expression(psi)

name <- paste(dirbase,"/",pname,"/nEX=",numEX,"/freq=",TLbase,"ps","/pmf/pmf_ADv",sep="")

cat(name,'\n')

#tiff(file.name,width=size.width,height=size.height)
#tiff(file.name,width=390,height=320)
width.size <- 390 * size
height.size <- 320 * size
file.name <- paste(name.out,'_',width.size,'-',height.size,'.tiff',sep='')
tiff(file.name,width=width.size,height=height.size)
fel3(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")
dev.off()
