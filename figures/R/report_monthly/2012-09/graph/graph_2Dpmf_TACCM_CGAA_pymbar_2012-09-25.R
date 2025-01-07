
source("~/calspa/GOLMAA/SH3/R/filled.contour.mod.R")
source("~/calspa/GOLMAA/SH3/R/fel3.R")
source("~/Rspa/plRamachad3.R")
source("~/Rspa/OutTiff.R")

dir0 <- "/home/yamamori/calspa/TACCM_CGAAREMD/AD"
dirbase <- paste(dir0,"/e_CG-FG_NH_2012-07-27",sep="")

nTZs<-1

name.title <- ""

title=name.title

label.x=expression(phi)
label.y=expression(psi)

name <- paste(dirbase,"/tau=",tau[1],"/mZ=",mZ[1],"/TZ=",TZs[1],"/",pname,"/nEX=",numEX,"/fq=",TLbase,"/pmf/pmf_pymbar_TAA=",TAA,"_TCG=",TCG,"_TZ=",TZs[1],"_KZAAo=",KZAAo,"_KZCGo=",KZCGo,"_",AACG,sep="")

cat(name,'\n')

file.name <- paste(name.out,"_",numuene,".tiff",sep='')

cat(file.name,'\n')
tiff(file.name,width=390,height=320)
fel3(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")
dev.off()

