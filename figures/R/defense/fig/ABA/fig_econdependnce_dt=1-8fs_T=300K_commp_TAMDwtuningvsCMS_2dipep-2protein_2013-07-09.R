source("~/Rspa/ini.R")

source("~/Rspa/set_enecon.sh")
source("~/Rspa/set_res.R")
source("~/Rspa/plmGeneral_wsrange.R")

leg.pos="topright"
is.leg <- 0

temp <- "300"

xrange <- c(0,10,10)
yrange <- c(-7,0,7)

name.res <- c("ASN","GLN")
name.pro <- c("BPTI","HP35")

num.res <- length(name.res)
num.pro <- length(name.pro)

numini=10

dirbaseABApept <- "/misc/ivy3/yamamori/ajisai/TAMD_econ/TAMD_econ/dipep/"
dirbaseABAprot="/misc/ivy3/yamamori/ajisai/TAMD/TAMD/mod_parm/proteins/"
dirbaseABAprot2="/misc/ivy3/yamamori/ajisai/TAMD_econ/TAMD_econ/proteins/"

file.name=paste("~/defense/eps/ABA/fig_econdependnce_dt=1-8fs_T=300K_commp_TAMDwtuningvsCMS_2dipep-2protein_2013-07-09.eps",sep='')
postscript(file.name,width=6,height=4.8,horizontal=FALSE,onefile=FALSE,paper="special")

hutosa <- 2

sen <- NULL
iro <- NULL
id.ys <- NULL
ids.ys <- NULL
tenshu <- NULL
senshu <- NULL

for ( k in 1:2  ) {
  sen[k]<-c(1)
  n<-1+k
  id.ys[k]<-c(2)
  ids.ys[k]<-c(3)
  tenshu[k]<-c(20)
}

iro[1]<-"red"
iro[2]<-"red"

senshu[1]<-1
senshu[2]<-2

tune <- c(0.1)
ntune <- length(tune)

par(mfrow=c(2,2))
par(oma=c(7.5,4.0,2.0,2.0))

file.names<-NULL
for (i in 1:num.res) {

  par(mar=c(0.0,0.0,0.0,0.0))

  file.names[1]=paste(dirbaseABApept,name.res[i],"/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_2012-04-17_pc6/dtune=",tune[1],"/",name.res[i],"Dv_econ_",temp[1],"_100ps_1-",numini,".econ.rmsd.av",sep='')
  file.names[2]=paste(dirbaseABApept,name.res[i],"/econ_wdi_NVE_AMBER_TermOn_afclustwoH_2012-04-17_pc6/",name.res[i],"Dv_econ_",temp[1],"_100ps_1-",numini,".econ.rmsd.av",sep='')
  
  plmGeneralwsrange(data.names=file.names,
                    sd.names=file.names,
                    id.ys=id.ys,
                    ids.ys=ids.ys,
                    label.size=0.5,axis.size=2.0,
                    iro=iro,axis.ft="F",is.header="T",
                    sdiro=iro,
                    xrange=xrange,yrange=yrange,
                    sdyrange=yrange,
                    is.sen=sen,width=10.0,
                    warrow="T")
  box(lwd=2.0)

  lines(c(1,6),c(-1.0,-1.0),lwd=2,lty="dashed",col="gray")
  
  if ( i==1 ) {
    yrange.axis <- c(-7,0,7)
    axis(2,yaxp=yrange.axis,lwd=2.0)
  }
}

file.names<-NULL
for (i in 1:num.pro) {

  par(mar=c(0.0,0.0,0.0,0.0))

  file.names[1]=paste(dirbaseABAprot,name.pro[i],"v","/econ_wdi_NVE_AMBER_TermOn_wtune_2012-03-26_pc6/","tune=",tune,"/",name.pro[i],"v_econ_",temp[1],"_100ps_1-",numini,".econ.rmsd.av",sep='')

  file.names[2]=paste(dirbaseABAprot2,name.pro[i],"/econ_wdi_NVE_AMBER_TermOn_afclust_2012-04-04_pc6/",name.pro[i],"v_econ_",temp[1],"_100ps_1-",numini,".econ.rmsd.av",sep='')

  plmGeneralwsrange(data.names=file.names,
                    sd.names=file.names,
                    id.ys=id.ys,
                    ids.ys=ids.ys,
                    label.size=0.5,axis.size=2.0,
                    iro=iro,axis.ft="F",is.header="T",
                    sdiro=iro,
                    xrange=xrange,yrange=yrange,
                    sdyrange=yrange,
                    is.sen=sen,width=10.0,
                    warrow="T")
  box(lwd=2.0)

  lines(c(1,6),c(-1.0,-1.0),lwd=2,lty="dashed",col="gray")

  if ( i==1 )   {
    xrange.axis <- c(0,9,9)
    axis(1,xaxp=xrange.axis,lwd=2.0)
  }
  if ( i==2 )   {
    xrange.axis <- c(0,10,10)
    axis(1,xaxp=xrange.axis,lwd=2.0)
  }
  
  if ( i==1 ) {
    yrange.axis <- c(-7,0,7)
    axis(2,yaxp=yrange.axis,lwd=2.0)
  }
}
