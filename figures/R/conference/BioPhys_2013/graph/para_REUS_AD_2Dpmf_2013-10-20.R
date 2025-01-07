
pname_REUS <-"UmbAt_Nbin_4x4_K_0.5"

numEX_REUS <-"500"

TLbase_REUS <-"10"

ff_REUS <-"ff99SB"

dir0 <- "~/calspa/refcalc/REUS/AD"
dirbase <- paste(dir0,"/s_REUSVAC_2013-08-30_",ff_REUS,sep="")
name <- paste(dirbase,"/",pname_REUS,"/nEX_",numEX_REUS,"/freq_",TLbase_REUS,"ps","/pmf/pmf_UmbMD_vac_nbx=40_nby=40.txt_2",sep="")

fact.x <- 180/pi
fact.y <- 180/pi
#fact.p <- 1.0/pi

dirout <- "/home/yamamori/gakkai/BioPhys_2013"

title <- NULL

label.x <- expression(paste(""))
label.y <- expression(paste(""))

name.out <- paste(dirout,"/eps/","fig_REUS_AD_2Dpmf_2013-10-20",sep='')
   
level <- seq(0.0,30.0,2.0)

xrange <- c(-180,180,6)
xrange.axis <- c(-180,180,6)
yrange <- c(-180,180,6)
yrange.axis <- c(-180,180,6)

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

