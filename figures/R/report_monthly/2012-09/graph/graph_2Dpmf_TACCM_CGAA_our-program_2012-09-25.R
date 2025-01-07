
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

name <- paste(dirbase,"/tau=",tau[1],"/mZ=",mZ[1],"/TZ=",TZs[1],"/",pname,"/nEX=",numEX,"/fq=",TLbase,"/pmf/pmf_TAA=",TAA,"_TCG_",TCG,"_TZ_",TZs[1],"_",width,"_",KZAAo,"_",KZCGo,"_",AACG,"_",numuene,"_2012-08-21",sep="")

cat(name,'\n')

file.name <- paste(name.out,"_",numuene,".tiff",sep='')

cat(file.name,'\n')
#tiff(file.name,width=390,height=320)
width.size <- 390 * size
height.size <- 320 * size
file.name <- paste(name.out,'_',width.size,'-',height.size,'.tiff',sep='')
tiff(file.name,width=width.size,height=height.size)
fel3(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")
dev.off()

