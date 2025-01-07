
source("ini.R")
source("set_spe.R")
source("set_res.R")
source("~/Rspa/plmGeneral_wsrange.R")

leg.pos="top"
is.leg <- 0

name.res <- c( "ASN", "GLN", "SER", "LYS" )
num.res <- length(name.res)

temp <- "300"

para.dihed.V <- "0.3"
num.para.dihed.V <- length(para.dihed.V)

xrange <- c(0.0,1.1,10)
xrange.axis <- c(0.1,1.0,9)

yrange <- c(-5,0.1,6)
yrange.axis <- c(-4.5,0.0,9)

num.para.dihed.V <- length(para.dihed.V)

par(mfrow=c(2,2))
par(oma=c(7.5,4.0,6.0,2.0))
par(mar=c(0.0,0.0,0.0,0.0))

dirbase.name <-paste("~/calspa/TAMD/mod_parm/dipep",sep='')
dirbase2.name <- paste("econ_wdi_NVE_AMBER_TermOn_wtune_2012-02-29_pc6",sep='')

label.names<-NULL
file.names<-NULL
title<-NULL

k <- 1
for ( i in 1:num.res  ) {
  j<-1

  cat(i,name.res[i],"\n")

  file.name <- paste(dirbase.name,"/",name.res[i],"/",dirbase2.name,"/",name.res[i],"Dv_econ_change_at_dt=5fs_300_100ps_1-10.econ.rmsd.av",sep='')

  id.xs<-NULL
  id.ys<-NULL
  ids.ys<-NULL
  tenshu<-NULL
  ave.x<-NULL
  sen<-NULL
  file.names <- NULL

  sen[1]<-c(1)                  
  iro[1] <- "red"

  id.xs[1] <- 1                 
  id.ys[1] <- 2
  ids.ys[1]<-c(3)
  tenshu[1]<-c(20)
  ave.x[1] <- 1
  file.names[1] <- file.name
  
  k <- k+1
  cat("k=",k,"\n")
  cat("id.xs=",id.xs,"\n")
  cat("id.ys=",id.ys,"\n")
  
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
  lines(c(0.1,1.0),c(-2.0,-2.0),lwd=2,lty="dashed",col="gray")
  
  if ( i == 3 || i == 4 ) {
    axis(1,xaxp=xrange.axis,lwd=2.0)
  }
  if ( i == 1 || i == 3 ) {
    axis(2,yaxp=yrange.axis,lwd=2.0)
  }
}

name.out=paste("~/seminars/GS/2012-02-29/fig4.econdependnce_atdt=5fs_tune_0.1-1.0_T=300K_commp_TAMDvsCMS_4dipep_2012-02-29",sep="")
OutTiff(name.out,w_val=600,h_val=500)
