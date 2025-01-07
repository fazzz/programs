
source("~/calspa/GOLMAA/SH3/R/filled.contour.mod.R")
source("~/calspa/GOLMAA/SH3/R/fel3.R")
source("~/Rspa/plRamachad3.R")

T    <-  c( "300"   ) 
T2    <- c(  300    ) 
nT   <- length(T)

TB    <-  c( "500"  )
nTB   <- length(TB)

KZ <- c( "300.0" )

mZ <- c( "50000.00" )

nmZ <- length(mZ)
nKZ <- length(KZ)

TAA <- "300"
TCG <- "370"
TZ <- "400"

p_g <- "0.50_3"
nEXT <- "1000"
width <- "0.1"

flag="CG"

level <- seq(-20,0,2)

label.x <- expression(paste(phi))
label.y <- expression(paste(psi))

for ( i in 1:nTB ) {
  for ( j in 1:nKZ ) {
    for ( k in 1:nmZ ) {
      name.title <- NULL

      title=name.title

      name <- paste("/home/yamamori/calspa/TACCM_CGAAREMD-s-ivy3/FiveAtomSysd/e_CG-AA_NH_2012-04-06/tau=1.0/mZ=50000.00/TZ=1000/ref/pmf/pmf_TAA=300_TCG_370_TZ_1000_0.1_",flag,sep="")

      cat(name,'\n')

      name.out <- paste("~/summary/CG-FG_TACCM_REMD/R/fig2_TACCM_CG_ref_2Dpmf_FAYSYS_2012-05-24",sep="")

      cat(name,'\n')
      
      fel3(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT")
      OutTiff(name.out,w_val=390,h_val=300)
    }
  }
}
