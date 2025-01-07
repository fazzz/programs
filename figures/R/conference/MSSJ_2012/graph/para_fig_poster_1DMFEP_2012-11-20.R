source("~/Rspa/plmGeneral_wsrange.R")

name.title <- NULL

name.out <- "~/gakkai/MSSJ_2012/tiff/fig_poster_1DMFEP_2012-11-20"

label.x=" "
label.y="pmf"

xrange <- c(0,1.0,10)
xrange.axis <- xrange
yrange <- c(0,6,6)
yrange.axis <- yrange

fact.x <- 1
fact.y <- 1
ave.x <- 1

file.name <- paste(name.out,".tiff",sep='')
cat(file.name,'\n')
tiff(file.name,width=300,height=200)

par(mfrow=c(1,1))
par(mar=c(0.0,0.0,0.0,0.0))
par(oma=c(3.0,6.0,2.0,2.0))

name <- "/home/yamamori/calspa/MFEP/AD/pmf2MGaussian_2012-10-10_bo_e_CG-FG_NH_2012-07-27/tau=1.0/mZ=100.00/TZ=700/SB_KZMAX=5000_NR4_woeljd0.001/nEX=1000/fq=10/path2_TAA=300_TCG=300_0_5000_CG_K=15_@-1.4,1.1-@1.1,-0.8.txt"

id.xs <- c(1,1)
id.ys <- c(2,2)
ids.ys <- c(3,3)
iro <- c(1,2)
senshu <- c(1,1)
tenshu <- c(3,3)
hutosa <- c(1,1)
is.leg <- 0

plmGeneralwsrange(data.names=name,
                  sd.names=name,
                  id.ys=id.ys,
                  ids.ys=ids.ys,
                  is.sen=rep(0,20),
                  label.size=0.5,axis.size=2.0,
                  iro=iro,axis.ft="F",is.header="T",
                  sdiro=iro,
                  xrange=xrange,yrange=yrange,
                  sdyrange=yrange,
                  warrow="T")
box(lwd=2.0)

axis(2,yaxp=yrange.axis,lwd=2.0,cex.axis=1.5)
mtext(outer=T,label.y,side=2,line=4.0,cex=1.5)

axis(1,xaxp=xrange.axis,lwd=2.0,cex.axis=1.5)
mtext(outer=T,label.x,side=1,line=4.0,cex=1.5)

