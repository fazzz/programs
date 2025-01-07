source("ini.R")
source("set_res.R")
source("plmGeneral_delta.R")

par(las=0)
xrange <- c(-5,5,10)
yrange <- c(0,0.5,5)
xrange3 <- c(0,100000,10)
yrange3 <- c(0,1.5,3)
xrangep <- c(-1.5,2,7)

is.leg <- 0

dirbase.name <-"/home/yamamori/work/programs/GIC/"
file.name1=paste(dirbase.name,'x.hist',sep='')
file.name2=paste(dirbase.name,'e.hist',sep='')
file.name3=paste(dirbase.name,'temp_mb.txt',sep='')
file.namep=paste(dirbase.name,'p.hist',sep='')
file.namef=paste(dirbase.name,'f.hist',sep='')

label.x1 <- "x"
label.x2 <- expression(eta)
label.x3 <- "step (/1000)"
label.xp <- "p"
label.xf <- "f"

label.y1 <- "f(x)"
label.y2 <- expression(paste("f(",eta,")"))
label.y3 <- "kT"
label.yp <- "f(p)"
label.yf <- "f(f)"

leg.pos <- "topright"

par(mfrow=c(3,2))
par(mar=c(0,5.0,6.0,2.0))
par(oma=c(7.5,2.0,2.0,2.0))

plmGeneralwsrange(data.names=file.name1,
                  id.xs=c(1),
                  id.ys=c(2),
                  iro=c(3),
                  xrange=xrange,
                  yrange=yrange,label.size=1.0,
                  is.sen=c(1),axis.ft="F",is.grid=1)
grid(col="lightgray",lty="dotted")
box(lwd=2.5)
axis(1,xaxp=xrange,lwd=2.5)
axis(2,yaxp=yrange,lwd=2.5)
mtext(outer=T,label.x1,side=1,line=3,cex=2.0)
mtext(outer=T,label.y1,side=2,line=3,cex=2.0)

plmGeneralwsrange(data.names=file.namep,
                  id.xs=c(1),
                  id.ys=c(2),
                  iro=c(3),
                  xrange=xrangep,
                  yrange=yrange,label.size=1.0,
                  is.sen=c(1),axis.ft="F",is.grid=1)
grid(col="lightgray",lty="dotted")
box(lwd=2.5)
axis(1,xaxp=xrange,lwd=2.5)
axis(2,yaxp=yrange,lwd=2.5)
mtext(outer=T,label.xp,side=1,line=3,cex=2.0)
mtext(outer=T,label.yp,side=2,line=3,cex=2.0)

plmGeneralwsrange(data.names=file.name2,
                  id.xs=c(1),
                  id.ys=c(2),
                  iro=c(3),
                  xrange=xrange,
                  yrange=yrange,label.size=1.0,
                  is.sen=c(1),axis.ft="F",is.grid=1)
grid(col="lightgray",lty="dotted")
box(lwd=2.5)
axis(1,xaxp=xrange,lwd=2.5)
axis(2,yaxp=yrange,lwd=2.5)
mtext(outer=T,label.x2,side=1,line=3,cex=2.0)
mtext(outer=T,label.y2,side=2,line=3,cex=2.0)

plmGeneralwsrange(data.names=file.namef,
                  id.xs=c(1),
                  id.ys=c(2),
                  iro=c(3),
                  xrange=xrangep,
                  yrange=yrange,label.size=1.0,
                  is.sen=c(1),axis.ft="F",is.grid=1)
grid(col="lightgray",lty="dotted")
box(lwd=2.5)
axis(1,xaxp=xrange,lwd=2.5)
axis(2,yaxp=yrange,lwd=2.5)
mtext(outer=T,label.xf,side=1,line=3,cex=2.0)
mtext(outer=T,label.yf,side=2,line=3,cex=2.0)

plmGeneralwsrange(data.names=file.name3,
                  id.xs=c(1),
                  id.ys=c(2),
                  iro=c(4),
                  xrange=xrange3,
                  yrange=yrange3,label.size=1.0,
                  is.sen=c(1),axis.ft="F",is.grid=1)
grid(col="lightgray",lty="dotted")
box(lwd=2.5)
axis(1,xaxp=xrange3,lwd=2.5)
axis(2,yaxp=yrange3,lwd=2.5)
mtext(outer=T,label.x3,side=1,line=6,cex=2.0)
mtext(outer=T,label.y3,side=2,line=3,cex=2.0)


name.out <- paste(dirbase.name,"distribution_x_eta",sep='')
OutTiff(name.out,w_val=900,h_val=1000)

