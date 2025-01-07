
name <- "~/calspa/refcalc/REUS/MetEnk/s_REUSVAC_2013-08-30_ff99SB/UmbAt_Nbin_8_K_0.4_2/nEX_1000/freq_10ps/pmf/pmf_MetEnkv"

fact.x <- 1
fact.y <- 1

dirout <- "/home/yamamori/gakkai/BioPhys_2013"

title <- NULL

label.x <- expression(paste(""))
label.y <- expression(paste(""))

name.out <- paste(dirout,"/eps/","fig_REUS_MetEnk_2Dpmf_2013-10-20",sep='')
   
level <- seq(0,20,2)

xrange <- c(0,0.3,3)       
xrange.axis <- c(0,0.3,3)  
yrange <- c(0,0.3,3)       
yrange.axis <- c(0.1,0.3,2)

file.name <- paste(name.out,'.eps',sep='')
postscript(file.name,width=4.0,height=2.5512,horizontal=FALSE,onefile=FALSE,paper="special")

par(oma = c(0,0,0,0) )
source("~/defense/fig/FillConMap.R")
source("~/defense/fig/FillConBar.R")
source("~/defense/fig/fel_FillConMap.R")
source("~/defense/fig/fel_FillConBar.R")

nf <- layout(matrix(c(1,2),1,2,byrow=TRUE),c(6.47,3.07),c(1))

par(cex.axis=1.2)
par(cex.lab=1.2)

cat(name)

label.x <- ""
label.y <- ""

felwFillConMap(name,label.x=label.x,label.y=label.y,level=level,
              title=title,kt="F/kBT",
              xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,
              mai=c(0.5,0.5,0.1,0.0))

felwFillConBar(name,label.x=label.x,label.y=label.y,level=level,
                   mai=c(1.348,0.75,0.1,0.75))

