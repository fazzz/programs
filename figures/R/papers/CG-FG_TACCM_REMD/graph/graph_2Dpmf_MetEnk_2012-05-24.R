 
source("~/calspa/GOLMAA/SH3/R/filled.contour.mod.R")
source("~/calspa/GOLMAA/SH3/R/fel3.R")
source("~/Rspa/plRamachad3.R")

dir0 <- "/home/yamamori/calspa/TACCM_CGAAREMD/MetEnk/"
dirbase <- paste(dir0,"/e_CG-FG_NH_2012-05-14/",state,sep="")

title=name.title

label.x="d1"
label.y="d2"

for ( i in 1:ntau ) {
  name <- paste(dirbase,"/tau=",tau[i],"/anl/pmf_MetEnk_AA_T=",TAA,"_",type,"_",width,".txt",sep="")

  cat(name,'\n')

  file.name <- paste(name.out,'.tiff',sep='')
  tiff(file.name,width=w_val,height=h_val)
  fel3(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")
  dev.off()
}
