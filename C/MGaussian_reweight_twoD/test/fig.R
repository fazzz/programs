
width <- "0.3"

fact.x=1
fact.y=1

level <- seq(0,10,0.5)

name.title <- NULL

MyColor <- function(n,alpha=1)
{			
    if ((n <- as.integer(n[1L])) > 0) {
        j <- n%/%3		
        k <- n%/%3		
        i <- n - j - k		
        c(if (i > 0) hsv(h = seq.int(from = 40/60, to = 25/60,
            length.out = i), alpha = alpha), if (j > 0) hsv(h = seq.int(from = 23/60,
            to = 11/60, length.out = j), alpha = alpha), if (k > 
            0) hsv(h = seq.int(from = 8/60, to = 0/60, length.out = k-1),
            alpha = alpha, s = seq.int(from = 1, to = 0.9, length.out = k-1),
            v = 1),hsv(0,0,1))
    }			
    else character(0L)	
}			

file.name <- "output/pmf12x12.txt.eps"
postscript(file.name,width=4.0,height=3.6,horizontal=FALSE,onefile=FALSE,paper="special")

minx<--1.0*pi
maxx<-pi
miny<--1.0*pi
maxy<-pi

xrange <- c(minx,maxx,6)
xrange.axis <- xrange
yrange <- c(miny,maxy,6)
yrange.axis <- yrange

par(mar = c(0.0,0.0,0.0,0.0) )
par(oma = c(0,0,0,0) )

source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange.R")
source("~/papers/CG-FG_TACCM_REMD/FillConBar.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap_wrange.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConBar_wmai.R")

nf <- layout(matrix(c(1,2),1,2,byrow=TRUE),c(3,1),1)

par(cex.axis=1.0)
par(cex.lab=1.0)

label.x <- expression(paste(phi,"(radian"))
label.y <- expression(paste(psi,"(radian"))

title=""
	
name <- "output/pmf12x12.txt"
cat(name)

felwFillConMapwrange(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange)

felwFillConBarwmai(name,label.x=label.x,label.y=label.y,level=level,mai=c(0.8,0.5,0.5,0.2))

