 
source("~/calspa/GOLMAA/SH3/R/filled.contour.mod.R")
source("~/calspa/GOLMAA/SH3/R/fel3.R")
source("~/Rspa/plRamachad3.R")

dir0 <- "/home/yamamori/calspa/TACCM_CGAAREMD/FiveAtomSysf"
dirbase <- paste(dir0,"/e_CG-FG_NH_2012-08-14_2",sep="")

nTZs<-1

title=name.title

label.x <- expression(paste(theta,1,sep=""))
label.y <- expression(paste(theta,2,sep=""))

name <- paste(dirbase,"/tau=",tau,"/mZ=",mZ,"/TZ=",TZs,"/",pname,"/nEX=",numEX,"/fq=",fq,"/pmf/pmf_TAA=",TAA,"_TCG_",TCG,"_TZ_",TZs,"_",width,"_",KZAAo,"_",KZCGo,"_",AACG,"_0",sep="")

cat(name,'\n')

file.name <- paste(name.out,'.tiff',sep='')
tiff(file.name,width=tiffwidth,height=tiffheight)
fel3(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")
dev.off()
