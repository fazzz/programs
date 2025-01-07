proname <- c( "CaM_apo_1cfd_mod" , "CaM_holo_1cll" )
nproname <- length(proname)

state <- c( "apo", "holo" )
nstate <- length(state)

T <- c( "450" , "500" )
nT <- 2

tau <- c( "1.0" , "1.0" )

cutoff <- c( "5.0", "6.5" )

ep <- c("0.50", "0.45")

name.title <- paste("/home/yamamori/Report/2012-12//tiff/fig_Q_single_well_CA_Go_Calmodulin_2013-01-14",sep="")

setwd("~/Rspa")

source("~/Rspa/ini.R")
source("~/Rspa/plmGeneral_wsrange_wdrange_wdiffx.R")

leg.pos="topright"
is.leg <- 0

xrange <- c(0,100000000,4)
xrange.axis <- c(10000000,100000000,3)

label.x <- expression(paste("time step "))

file.name <- paste(name.title,".tiff",sep="")
tiff(file.name,width=800,height=400)

par(mfrow=c(1,1))
par(mar=c(0.0,0.0,0.0,0.0))
par(oma=c(7.5,6.0,2.0,2.0))

filex.names<-NULL
filey.names<-NULL
#yrange <- c(-0.1,1.1,10)
yrange <- c(-0.01,1.01,10)
yrange.axis <- c(0,1,10)
label.y <- expression(paste("Q"))

id.ys<-NULL

hutosa <- 1

dir <- "~/calspa/CaM/Go-like_model/SB/e_vT_vpara_GOLMCA_JMB2001_NH_2012-07-19/"

sen<-1
senshu <- 1
iro<-2
id.ys<-c(1)

for (i in 1:1) {
    dirbase.name <- paste(dir,"/",state[i],"/tau=",tau,"/ep=",ep[i],"/coff=",cutoff[i],"/",sep="")
        
    filex.names <- paste(dirbase.name,proname[i],"_equ_T=",T[i],".out",sep='')
    filey.names <- paste(dirbase.name,"anl/",proname[i],"_equ_T=",T[i],".qca",sep='')
    ave.x <- 1

    plmGeneralwsrangewdrangewdiffx(datax.names=filex.names,
                                   datay.names=filey.names,
                                   sd.names=file.names,
                                   id.ys=id.ys,
                                   label.size=0.5,axis.size=2.0,
                                   iro=iro,axis.ft="F",is.header="F",
                                   sdiro=iro,
                                   xrange=xrange,yrange=yrange,
                                   sdyrange=yrange,
                                   is.sen=sen,width=10.0,
                                   warrow="F",ave.x=ave.x)
   box(lwd=2.0)
          
   axis(1,xaxp=xrange.axis,lwd=2.0,cex=2.0)
   mtext(outer=T,label.x,side=1,line=5.0,cex=2.0)
          
   if ( i == 1 ) {
       axis(2,yaxp=yrange.axis,lwd=2.0,cex=2.0)
       mtext(label.y,side=2,line=3.0,cex=2.0)
   }
}

dev.off()

