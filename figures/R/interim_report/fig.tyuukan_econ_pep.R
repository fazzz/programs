source("~/Rspa/ini.R")
source("~/Rspa/set_enecon.sh")
source("~/Rspa/plmGeneral_wsrange.R")

is.leg <- 0

xrange <- c(0,6,6)
yrange <- c(-6,0,6)

temp <- 300
dt <- c("1","2","3","4","5","6","7","8","9")
num.temp <- length(temp)
num.dt <- length(dt)

prot.name <- c("GLY20","villin head piece subdomain")
#prot.name <- c("GLY20","villin head piece subdomain")
num.prot <- length(prot.name)

par(mfrow=c(1,num.prot))
par(mar=c(6.0,0.0,0.0,0.0))
par(oma=c(7.5,5.0,2.0,2.0))

for ( k in 1:4  ) {
  sen[k]<-c(3)
  n<-1+k
  iro[k]<-c(n)
  id.ys[k]<-c(2)
  tenshu[k]<-c(20)
}

file.names <- NULL
label.names[1]<-"pc"
label.names[2]<-"vV"
label.names[3]<-"CMD"

for ( i in 1:num.prot ) {
  if (i==1) {
    file.names[1]<-"/home/yamamori/calspa/TAMD/decapeptides/GLY_20/econ_wdi_NVE_amber_100ps_termon_dt_2_2/GLY_20v_econ_300_100ps_1-10.econ.rmsd.av"
    file.names[2]<-"/home/yamamori/calspa/TAMD/decapeptides/GLY_20/econ_wdi_NVE_amber_100ps_termon_dt_2_ai_debug/GLY_20v_econ_300_100ps_1-10.econ.rmsd.av"
    file.names[3]<-"/home/yamamori/calspa/refcalc/G20/econ_wdi_NVE_100ps/GLY_20v_econ_300_100ps_1-10_2ntc.econ.rmsd.av"
  }
  else if (i==2) {
    file.names[1]<-"/home/yamamori/calspa/TAMD/HP35/econ_wdi_NVE_amber_100ps_termon_dt_2_2/HP35v_econ_300_100ps_1-10.econ.rmsd.av"    
    file.names[2]<-"/home/yamamori/calspa/TAMD/HP35/econ_wdi_NVE_amber_100ps_termon_dt_2_ai_debug/HP35v_econ_300_100ps_1-3.econ.rmsd.av"
    file.names[3]<-"/home/yamamori/calspa/refcalc/HP35/econ_wdi_NVE_100ps/HP35v_econ_300_100ps_1-10_2ntc.econ.rmsd.av"
  }

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
  text(3.0,-5.5,prot.name[i],cex=1.5)
  box(lwd=2.0)
  a <- read.table(file.names[3],header=TRUE)
  cat(a[2,2])
  lines(c(1,6),c(a[2,2],a[2,2]),lwd=2,lty="dashed",col="gray")
  
  axis(1,xaxp=c(0,5,5),lwd=2.0,cex.axis=1.5)
#  mtext(outer=T,label.x,side=1,line=5,cex=1.5)

  if ( i==1 ) {
    axis(2,yaxp=yrange,lwd=2.0,cex.axis=1.5)
#    mtext(outer=T,label.y,side=2,line=3.0,cex=1.5)
  }

}

name.out=paste("~/papers/TAMD/fig.tyuukan_econ_pep",sep="")
OutTiff(name.out,w_val=800,h_val=500)

