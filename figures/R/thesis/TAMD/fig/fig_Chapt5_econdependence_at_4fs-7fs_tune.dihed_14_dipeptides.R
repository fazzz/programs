
source("~/Rspa/ini.R")
source("~/Rspa/set_spe.R")
source("~/Rspa/set_res.R")
source("~/Rspa/plmGeneral_wsrange.R")

leg.pos="top"
is.leg <- 0

name.res <- c( "ASN", "GLN", "SER", "LYS", "ASP", "ARG", "CYS", "GLU", "ILE", "LEU", "LYS", "MET", "TRP", "THR", "VAL" )
num.res <- length(name.res)

temp <- "300"

hutosa <- 1

para.dihed.V <- "0.3"
num.para.dihed.V <- length(para.dihed.V)

xrange <- c(0.0,1.1,10)
xrange.axis <- c(0.1,1.0,9)

yrange <- c(-5,0.1,6)
yrange.axis <- c(-4.5,0.0,9)

num.para.dihed.V <- length(para.dihed.V)

#dt <- c("5","6","7")
dt <- c("7","6","5")
num.dt <- length(dt)

name.out=paste("~/thesis/TAMD/eps/fig_Chapt5_econdependence_at_4fs-7fs_tune.dihed_14_dipeptides.eps",sep="")
postscript(name.out,width=6.0,height=6.0,horizontal=FALSE,onefile=FALSE,paper="special")

par(mfrow=c(3,5))
par(oma=c(5.5,5.5,2.0,2.0))
par(mar=c(0.0,0.0,0.0,0.0))

dirbase.name <-paste("~/calspa/TAMD-s-ivy3/TAMD/mod_parm/dipep",sep='')
dirbase2.name <- paste("econ_wdi_NVE_AMBER_TermOn_wtune_2012-02-29_pc6",sep='')

label.names<-NULL
file.names<-NULL
title<-NULL

iro <- c("red","blue","green")

k <- 1
for ( i in 1:num.res  ) {
  j<-1
  cat(i,name.res[i],"\n")
  id.xs<-NULL
  id.ys<-NULL
  ids.ys<-NULL
  tenshu<-NULL
  ave.x<-NULL
  sen<-NULL
  file.names <- NULL
  

  for ( t in 1:num.dt  ) {
    file.names[t]<- paste(dirbase.name,"/",name.res[i],"/",dirbase2.name,"/",name.res[i],"Dv_econ_change_at_dt=",dt[t],"fs_300_100ps_1-10.econ.rmsd.av",sep='')

    sen[t]<-c(1)                  
    id.xs[t] <- 1                 
    id.ys[t] <- 2
    ids.ys[t]<-c(3)
    tenshu[t]<-c(20)
    ave.x[t] <- 1
  }
  
  
  plmGeneralwsrange(data.names=file.names,
                    sd.names=file.names,
                    id.xs=id.xs,
                    id.ys=id.ys,
                    ave.x=ave.x,
                    ids.ys=ids.ys,
                    label.size=0.5,axis.size=2.0,
                    iro=iro,axis.ft="F",is.header="T",
                    sdiro=iro,
                    xrange=xrange,yrange=yrange,
                    sdyrange=yrange,
                    is.sen=sen,
                    width=10.0,
                    warrow="T")
  box(lwd=2.0)
  text(0.5,-0.3,name.res[i])
  lines(c(0.1,1.0),c(-2.0,-2.0),lwd=2,lty="dashed",col="gray")
  
  if ( i > 10 ) {
    axis(1,xaxp=xrange.axis,lwd=2.0)
  }

  if ( i == 1 || i == 6 || i == 11  ) {              
    axis(2,yaxp=yrange.axis,lwd=2.0)
  }                                 
  
}

label.x <- expression(paste("time step (fs)"))
mtext(outer=T,label.x,side=1,line=3.0,cex=1.0)        
label.y <- expression(paste("energy (kcal/mol)"))
mtext(outer=T,label.y,side=2,line=3.0,cex=1.0)
