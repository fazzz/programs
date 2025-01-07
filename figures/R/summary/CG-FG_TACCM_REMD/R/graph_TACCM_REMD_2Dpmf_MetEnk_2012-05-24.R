 
source("~/calspa/GOLMAA/SH3/R/filled.contour.mod.R")
source("~/calspa/GOLMAA/SH3/R/fel3.R")
source("~/Rspa/plRamachad3.R")

dir0 <- "/home/yamamori/calspa/TACCM_CGAAREMD/MetEnk/"
dirbase <- paste(dir0,"/e_CG-FG_NH_2012-05-14/",state,sep="")

nTZs<-1

title=name.title

label.x="d1"
label.y="d2"

for ( i in 1:ntau ) {
  for ( j in 1:nmZ ) {
    for ( k in 1:nTZs ) {
      name <- paste(dirbase,"/tau=",tau[i],"/mZ=",mZ[j],"/ep=",ep,"/cutoff=",cutoff,"/TZ=",TZs[k],"/",pname,"/nEXT=",numEX,"/pmf/pmf_TAA=",TAA,"_TCG1_",TCG1,"_TCG2_",TCG2,"_TZ_",TZs[k],"_",width,"_",KZAAo,"_",KZCG1o,"_",KZCG2o,"_",AACG,sep="")

      cat(name,'\n')

      fel3(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")
    }
  }
}


