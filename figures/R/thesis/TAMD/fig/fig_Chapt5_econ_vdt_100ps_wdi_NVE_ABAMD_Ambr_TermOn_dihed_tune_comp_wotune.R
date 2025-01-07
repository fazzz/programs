source("~/Rspa/ini.R")

source("~/Rspa/set_enecon.sh")
source("~/Rspa/set_res.R")
source("~/Rspa/plmGeneral_wsrange.R")

leg.pos="topright"
is.leg <- 0

temp <- "300"

hutosa <- 1

xrange <- c(0,8,8)
xrange.axis <- c(0,7,7)
yrange <- c(-7,0,7)
yrange <- c(-7,0,7)

numini=10

name.proteins <- c("UBIQ", "BPTI", "HP35")
num.proteins <- length(name.proteins)
tune <- c(0.1)
ntune <- length(tune)

dirbase.name <- "/home/yamamori/calspa/TAMD-s-ivy3/TAMD/mod_parm/proteins/"
dirbase2.name <- "/home/yamamori/calspa/TAMD_econ-s-ivy3/TAMD_econ/proteins/"

name.out=paste("~/thesis/TAMD/eps/fig_Chapt5_econ_vdt_100ps_wdi_NVE_ABAMD_Ambr_TermOn_dihed_tune_comp_wotune.eps",sep='')
postscript(name.out,width=4.0,height=3.5,horizontal=FALSE,onefile=FALSE,paper="special")

id.ys <- NULL
ids.ys <- NULL

for ( k in 1:4  ) {
  sen[k]<-c(1)
  n<-1+k
  iro[k]<-c(n)
  id.ys[k]<-c(2)
  ids.ys[k]<-c(3)
  tenshu[k]<-c(20)
}

iro[1]<-"red"
iro[2]<-"blue"
iro[3]<-"green"

par(mfrow=c(1,num.proteins))
par(oma=c(7.5,4.0,2.0,2.0))
par(mar=c(0.0,0.0,0.0,0.0))

file.names<-NULL
k2<-1
for (i in 1:num.proteins) {
  for (j in 1:ntune) {
    file.names[1]=paste(dirbase.name,name.proteins[i],"v","/econ_wdi_NVE_AMBER_TermOn_wtune_2012-03-26_pc6/","tune=",tune,"/",name.proteins[i],"v_econ_",temp[1],"_100ps_1-",numini,".econ.rmsd.av",sep='')

    file.names[2]=paste(dirbase2.name,name.proteins[i],"/econ_wdi_NVE_AMBER_TermOn_afclust_2012-04-04_pc6/",name.proteins[i],"v_econ_",temp[1],"_100ps_1-",numini,".econ.rmsd.av",sep='')

    file.names[3]=paste("/home/yamamori/calspa/refcalc-s-ivy3/refcalc/",name.proteins[i],"/econ_wdi_NVE_100ps/",name.proteins[i],"v_econ_",temp[1],"_100ps_1-10_2ntc",".econ.rmsd.av",sep='')

    senshu[1]<-1
    senshu[2]<-1
  
    plmGeneralwsrange(data.names=file.names,
                      sd.names=file.names,
                                        #                      label.names=label.names,
                      id.ys=id.ys,
                      ids.ys=ids.ys,
                      label.size=0.5,axis.size=2.0,
                      iro=iro,axis.ft="F",is.header="T",
                      sdiro=iro,
                      xrange=xrange,yrange=yrange,
                      sdyrange=yrange,
                      is.sen=sen,width=10.0,
                      ,warrow="T")
    text(5,-6,name.proteins[i])
    box(lwd=2.0)
  
#  a <- read.table(file.names[3],header=TRUE)
#  cat(a[2,2])
    lines(c(1,6),c(-1,-1),lwd=2,lty="dashed",col="gray")
    box(lwd=2.0)

    axis(1,xaxp=xrange.axis,lwd=2.0)
    mtext(outer=T,label.x,side=1,line=5.0,cex=0.8)
  
    if (i==1 ) {
      axis(2,yaxp=yrange,lwd=2.0)
      mtext(outer=T,label.y,side=2,line=2.5,cex=0.8)
    }
  }
}

