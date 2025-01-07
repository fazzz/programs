
fact.x <- 180/pi
fact.y <- 180/pi

dirout <- "/home/yamamori/gakkai/News_Letter_2014"

title <- NULL

label.x <- expression(paste(""))
label.y <- expression(paste(""))

name.out <- paste(dirout,"/eps/","fig_AD_2Dpmf_2014-02-05",sep='')
   
level <- seq(0.0,30.0,2.0)

xrange <- c(-180,180,6)
xrange.axis <- c(-180,180,6)
yrange <- c(-180,180,6)
yrange.axis <- c(-180,180,6)

file.name <- paste(name.out,'.eps',sep='')
postscript(file.name,width=2.7,height=5.6,horizontal=FALSE,onefile=FALSE,paper="special")

par(oma = c(0,0,0,0) )
source("~/defense/fig/FillConMap.R")
source("~/defense/fig/FillConBar.R")
source("~/defense/fig/fel_FillConMap.R")
source("~/defense/fig/fel_FillConBar.R")

nf <- layout(matrix(c(1,6,2,6,3,6,4,6,5,6),5,2,byrow=TRUE),c(3,1.2),c(20,15,15,15,23))

par(cex.axis=1.2)
par(cex.lab=1.2)

label.x <- ""
label.y <- ""

name<-paste("/home/yamamori/data_paper/MuSTAR_MD/fig_3/pmf_AD_MuSTAR_MD",sep='')
cat(name)
felwFillConMap(name,label.x=label.x,label.y=label.y,level=level,
              title=title,kt="F/kBT",
              xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,
              xaflag="F",
              mai=c(0.0,0.5,0.3,0.1))
text(120,120,"(a)",cex=1.6)

name<-paste("/home/yamamori/data_paper/MuSTAR_MD/fig_3/pmf_AD_TAMD",sep='')
cat(name)
felwFillConMap(name,label.x=label.x,label.y=label.y,level=level,
              title=title,kt="F/kBT",
              xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,
              xaflag="F",
              mai=c(0.0,0.5,0.0,0.1))
text(120,120,"(b)",cex=1.6)

name<-paste("/home/yamamori/data_paper/MuSTAR_MD/fig_3/pmf_AD_REMD",sep='')
cat(name)
felwFillConMap(name,label.x=label.x,label.y=label.y,level=level,
              title=title,kt="F/kBT",
              xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,
              xaflag="F",
              mai=c(0.0,0.5,0.0,0.1))
text(120,120,"(c)",cex=1.6)

name<-paste("/home/yamamori/data_paper/MuSTAR_MD/fig_3/pmf_AD_REUS",sep='')
cat(name)
felwFillConMap(name,label.x=label.x,label.y=label.y,level=level,
              title=title,kt="F/kBT",
              xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,
              xaflag="F",
              mai=c(0.0,0.5,0.0,0.1))
text(120,120,"(d)",cex=1.6)

name<-paste("/home/yamamori/data_paper/MuSTAR_MD/fig_3/pmf_AD_CMD",sep='')							   
cat(name)										   
felwFillConMap(name,label.x=label.x,label.y=label.y,level=level,			   
              title=title,kt="F/kBT",							   
              xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis, 
              mai=c(0.5,0.5,0.0,0.1))							   
text(120,120,"(e)",cex=1.6)

felwFillConBar(name,label.x=label.x,label.y=label.y,level=level,
               mai=c(2.0,0.3,0.3,0.6))

