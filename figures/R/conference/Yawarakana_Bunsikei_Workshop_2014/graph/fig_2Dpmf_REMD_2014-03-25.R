
fact.x <- 180/pi
fact.y <- 180/pi

dirout <- "/home/yamamori/gakkai/Yawarakana_Bunsikei_Workshop_2014/"

title <- NULL

label.x <- expression(paste(""))
label.y <- expression(paste(""))

name.out <- paste(dirout,"/eps/","fig_AD_2Dpmf_REMD_2014-03-25",sep='')
   
level <- seq(0.0,30.0,2.0)

xrange <- c(-180,180,6)
xrange.axis <- c(-180,180,6)
yrange <- c(-180,180,6)
yrange.axis <- c(-180,180,6)

file.name <- paste(name.out,'.eps',sep='')
postscript(file.name,width=4.1496,height=2.6472,horizontal=FALSE,onefile=FALSE,paper="special")

par(oma = c(0,0,0,0) )
source("~/defense/fig/FillConMap.R")
source("~/defense/fig/FillConBar.R")
source("~/defense/fig/fel_FillConMap.R")
source("~/defense/fig/fel_FillConBar.R")

nf <- layout(matrix(c(1,2),1,2,byrow=TRUE),c(3,1.2),c(1))

par(cex.axis=1.2)
par(cex.lab=1.2)

label.x <- ""
label.y <- ""

name<-paste("/home/yamamori/data_paper/MuSTAR_MD/fig_3/pmf_AD_REMD",sep='')
cat(name)
felwFillConMap(name,label.x=label.x,label.y=label.y,level=level,
              title=title,kt="F/kBT",
              xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,
              mai=c(0.5,0.5,0.1,0.0))

felwFillConBar(name,label.x=label.x,label.y=label.y,level=level,
               mai=c(1.348,0.75,0.1,0.75))

