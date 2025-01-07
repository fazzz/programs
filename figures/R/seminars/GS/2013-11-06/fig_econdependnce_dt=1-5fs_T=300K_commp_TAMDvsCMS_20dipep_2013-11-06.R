source("~/Rspa/ini.R")

source("~/Rspa/set_enecon.sh")
source("~/Rspa/set_res.R")
source("~/Rspa/plmGeneral_wsrange.R")

leg.pos="topright"
is.leg <- 0

temp <- "300"

xrange <- c(0,10,10)
yrange <- c(-7,0,7)

name.res <- c( "ASN", "GLN", "LEU","SER", "LYS","ILE","ALA", "CYS", "MET", "THR", "VAL", "ARG", "ASP", "GLY", "HIE","PHE","PRO", "THR","TYR","TRP" )

num.res <- length(name.res)

numini=10

dirbaseABA <- "/misc/ivy3/yamamori/ajisai/TAMD_econ/TAMD_econ/dipep/"
dirbaseCMD="/misc/ivy3/yamamori/ajisai/refcalc/refcalc/dipep/"

file.name=paste("~/seminars/GS/2013-11-06/fig_econdependnce_dt=1-5fs_T=300K_commp_TAMDvsCMS_20dipep_2013-11-06.eps",sep='')
postscript(file.name,width=8,height=6,horizontal=FALSE,onefile=FALSE,paper="special")

hutosa <- 1

sen <- NULL
iro <- NULL
id.ys <- NULL
ids.ys <- NULL
tenshu <- NULL
senshu <- NULL

for ( k in 1:2  ) {
  sen[k]<-c(1)
  n<-1+k
  iro[k]<-c(n)
  id.ys[k]<-c(2)
  ids.ys[k]<-c(3)
  tenshu[k]<-c(20)
  senshu[k]<-1
}

par(mfrow=c(4,5))
par(oma=c(7.5,4.0,2.0,2.0))

file.names<-NULL
for (i in 1:num.res) {

  par(mar=c(0.0,0.0,0.0,0.0))

  file.names[1]=paste(dirbaseABA,name.res[i],"/econ_wdi_NVE_AMBER_TermOn_afclustwoH_2012-04-17_pc6/",name.res[i],"Dv_econ_",temp[1],"_100ps_1-",numini,".econ.rmsd.av",sep='')
  file.names[2]=paste(dirbaseCMD,name.res[i],"/econ_wdi_NVE_100ps/",name.res[i],"Dv_econ_",temp[1],"_100ps_1-10_2ntc",".econ.rmsd.av",sep='')
  
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
  text(6,-6.0,name.res[i],cex=1.5)
  box(lwd=2.0)

  lines(c(1,6),c(-1.0,-1.0),lwd=2,lty="dashed",col="gray")

  if (i==16 || i==17 || i==18 || i==19 ) {
    xrange.axis <- c(0,9,9)
    axis(1,xaxp=xrange.axis,lwd=2.0)
  }

  if (i==20) {
    xrange.axis <- c(0,10,10)
    axis(1,xaxp=xrange.axis,lwd=2.0)
  }

  if (i==1 ) {
    yrange.axis <- c(-7,0,7)
    axis(2,yaxp=yrange.axis,lwd=2.0)
  }
  
  if (i==6 || i==11 || i==16 ) {
    yrange.axis <- c(-7,-1,6)
    axis(2,yaxp=yrange.axis,lwd=2.0)
  }
  
}
