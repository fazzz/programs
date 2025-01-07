source("~/Rspa/ini.R")
source("~/Rspa/set_enecon.sh")
source("~/Rspa/plmGeneral_wsrange.R")

is.leg <- 1

xrange <- c(0,6,6)
yrange <- c(-6,0,6)

temp <- 300
dt <- c("1","2","3","4","5","6","7","8","9")
num.temp <- length(temp)
num.dt <- length(dt)

prot.name <- c("ubiqutin")
num.prot <- length(prot.name)

par(mfrow=c(1,num.prot))
par(mar=c(4.0,0.0,0.0,0.0))
par(oma=c(5.0,5.0,2.0,2.0))

for ( k in 1:4  ) {
  sen[k]<-c(3)
  n<-2+k
  iro[k]<-c(n)
  id.ys[k]<-c(2)
  tenshu[k]<-c(20)
}

file.names <- NULL
label.names[1]<-"CMD"
label.names[2]<-"pc"
label.names[3]<-"vV"

for ( i in 1:num.prot ) {
  file.names[1]<-"/home/yamamori/calspa/refcalc/UBIQ/econ_wdi_NVE_100ps/UBIQv_econ_300_100ps_1-10_2ntc.econ.rmsd.av"
  file.names[2]<-"/home/yamamori/calspa/TAMD/UBIQ/econ_wdi_NVE_amber_100ps_termon_dt_2_2/UBIQv_econ_300_100ps_1-10.econ.rmsd.av"
  file.names[3]<-"/home/yamamori/calspa/TAMD/UBIQ/econ_wdi_NVE_amber_100ps_termon_dt_2_2_ai_debug/UBIQv_econ_300_100ps_1-10.econ.rmsd.av"

  plmGeneralwsrange(data.names=file.names,
                    sd.names=file.names,
                    label.names=label.names,
                    id.ys=id.ys,
                    ids.ys=ids.ys,
                    label.size=0.5,axis.size=2.0,
                    iro=iro,axis.ft="F",is.header="T",
                    sdiro=iro,
                    xrange=xrange,yrange=yrange,
                    sdyrange=yrange,
                    is.sen=sen,width=10.0,
                    ,warrow="T")
  text(3.0,-5.5,prot.name[i])
  box(lwd=2.0)
  a <- read.table(file.names[1],header=TRUE)
  cat(a[2,2])
  lines(c(1,6),c(a[2,2],a[2,2]),lwd=2,lty="dashed",col=iro[4])
  
  axis(1,xaxp=xrange,lwd=2.0)
  mtext(outer=T,label.x,side=1,line=5,cex=1.5)

  axis(2,yaxp=yrange,lwd=2.0)
  mtext(outer=T,label.y,side=2,line=3.0,cex=1.5)

}

name.out=paste("~/papers/TAMD/fig_abstract(enecon_UBIQ)",sep="")
OutTiff(name.out,w_val=500,h_val=500)

