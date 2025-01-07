 
source("~/calspa/GOLMAA/SH3/R/filled.contour.mod.R")
source("~/calspa/GOLMAA/SH3/R/fel3.R")
source("~/Rspa/plRamachad3.R")

dir0 <- "/home/yamamori/calspa/TACCM_CGAAREMD/AD"
dirbase <- paste(dir0,"/e_CG-FG_NH_2012-05-12",sep="")

nTZs<-1

title=name.title

label.x <- expression(paste(phi))
label.y <- expression(paste(psi))

for ( i in 1:ntau ) {
  for ( j in 1:nmZ ) {
    for ( k in 1:nTZs ) {
      name <- paste(dirbase,"/tau=",tau[i],"/mZ=",mZ[j],"/TZ=",TZs[k],"/",pname,"/nEX=",numEX,"/pmf/pmf_TAA=",TAA,"_TCG_",TCG,"_TZ_",TZs[k],"_",width,"_",KZAAo,"_",KZCGo,"_",AACG,sep="")

      cat(name,'\n')

      file.name <- paste(name.out,'.tiff',sep='')
      tiff(file.name,width=390,height=320)
      fel3(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")
      dev.off()
    }
  }
}
