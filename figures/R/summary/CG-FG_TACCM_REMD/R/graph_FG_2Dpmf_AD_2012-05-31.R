 
source("~/calspa/GOLMAA/SH3/R/filled.contour.mod.R")
source("~/calspa/GOLMAA/SH3/R/fel3.R")
source("~/Rspa/plRamachad3.R")

dir0 <- "/home/yamamori/calspa/TACCM_CGAAREMD/AD"
dirbase <- paste(dir0,"/e_FG_NH_2012-05-31",sep="")

nTZs<-1

title=name.title

label.x <- expression(paste(phi))
label.y <- expression(paste(psi))

for ( i in 1:ntau ) {
  name <- paste(dirbase,"/tau=",tau[i],"/anl/pmf_AD_AA_T=",TAA,"_",width,"_",TSns,"ns.txt",sep="")

  cat(name,'\n')

  file.name <- paste(name.out,'.tiff',sep='')
  tiff(file.name,width=390,height=320)
  fel3(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")
  dev.off()
}
