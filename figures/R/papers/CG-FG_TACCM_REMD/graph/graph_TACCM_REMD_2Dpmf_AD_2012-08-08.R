 
source("~/calspa/GOLMAA/SH3/R/filled.contour.mod.R")
source("~/calspa/GOLMAA/SH3/R/fel3.R")
source("~/Rspa/plRamachad3.R")

dir0 <- "/home/yamamori/calspa/TACCM_CGAAREMD/AD"
dirbase <- paste(dir0,"/e_CG-FG_NH_2012-07-27",sep="")

nTZs<-1

title=name.title

label.x <- expression(paste(phi))
label.y <- expression(paste(psi))

name <- paste(dirbase,"/tau=",tau,"/mZ=",mZ,"/TZ=",TZ,"/",pname,"/nEX=",numEX,"/fq=",fq,"/pmf/pmf_TAA=",TAA,"_TCG_",TCG,"_TZ_",TZ,"_",width,"_",KZAAo,"_",KZCGo,"_",AACG,"_4_2012-08-21",sep="")

cat(name,'\n')

title <- NULL

file.name <- paste(name.out,'.tiff',sep='')
tiff(file.name,width=390,height=320)
fel3(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")
dev.off()
