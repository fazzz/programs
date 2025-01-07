
source("ini.R")
source("set_spe.R")
source("set_res.R")
source("~/Rspa/plmGeneral_wsrange.R")

leg.pos="top"
is.leg <- 0

name.res <- c( "ASN", "GLN" )
num.res <- length(name.res)

temp <- "300"

para.dihed.V <- "0.3"
num.para.dihed.V <- length(para.dihed.V)

xrange <- c(0.09,1.01,10)
xrange.axis <- c(0.1,1.0,9)
#xrange <- c(-100,100,10)

yrange <- c(-6,1,7)
yrange.axis <- c(-6,1,7)
#yrange <- c(-100,100,10)

num.para.dihed.V <- length(para.dihed.V)

dt <- c("7","6","5")
#dt <- c("4")
num.dt <- length(dt)

par(mfrow=c(1,2))
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

    if ( t==1)
      iro[t]="red"
    if ( t==2)
      iro[t]="blue"
    if ( t==3)
      iro[t]="green"
      
    id.xs[t] <- 1                 
    id.ys[t] <- 2
    ids.ys[t]<- 3
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
#  lines(c(0.1,1.0),c(-2.0,-2.0),lwd=2,lty="dashed",col="gray")
  
  axis(1,xaxp=xrange.axis,lwd=2.0)
  
  if ( i == 1  ) {               
    axis(2,yaxp=yrange.axis,lwd=2.0) 
  }                                  
  
}

name.out=paste("~/gakkai/seibutu_kanto/fig2_econdependence_at_4fs-7fs_tune.dihed_ASN_GLN_dipeptides_2012-03-01",sep="")
OutTiff(name.out,w_val=600,h_val=500)
