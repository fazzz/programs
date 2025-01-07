source("~/Rspa/ini.R")

source("~/Rspa/set_enecon.sh")
source("~/Rspa/set_res.R")
source("~/Rspa/plmGeneral_wsrange.R")

leg.pos="topright"
is.leg <- 0

temp <- "300"

xrange <- c(0,9,9)
yrange <- c(-7,0.0,8)
yrange.axis <- c(-7,-1,6)

lambda <- "0.1"

name.res <- c( "ASN", "GLN", "LYS", "SER", "ALA", "ARG", "CYS", "ILE", "LEU", "MET", "THR", "VAL" )

numini=10

num.res <- length(name.res)

dirbase.name <- "/home/yamamori/calspa/TAMD/dipeptides/"
dirbase21.name <-paste("~/calspa/TAMD/mod_parm/dipep",sep='')
dirbase22.name <- paste("econ_wdi_NVE_AMBER_TermOn_wtune_2012-02-29_pc6",sep='')

sen<-NULL
id.ys<-NULL
ids.ys<-NULL
tenshu<-NULL

for ( k in 1:3  ) {
  sen[k]<-c(1)
  id.ys[k]<-c(2)
  ids.ys[k]<-c(3)
  tenshu[k]<-c(20)
}

iro[1] <- "red"
iro[2] <- "black"
iro[3] <- "blue"

par(mfrow=c(4,5))
par(oma=c(7.5,4.0,2.0,2.0))
par(mar=c(0.0,0.0,0.0,0.0))

file.names<-NULL
k2<-1
for (i in 1:num.res) {
  file.names[1]<- paste(dirbase21.name,"/",name.res[i],"/",dirbase22.name,"/tune=",lambda,"/",name.res[i],"Dv_econ_",temp[1],"_100ps_1-",numini,".econ.rmsd.av",sep='')
  file.names[2]<- paste(dirbase21.name,"/",name.res[i],"/",dirbase22.name,"/tune=1.0/",name.res[i],"Dv_econ_",temp[1],"_100ps_1-",numini,".econ.rmsd.av",sep='')

  file.names[3]=paste("/home/yamamori/calspa/refcalc/dipep/",name.res[i],"/econ_wdi_NVE_100ps/",name.res[i],"Dv_econ_",temp[1],"_100ps_1-10_2ntc",".econ.rmsd.av",sep='')

  senshu[1]<-1
  senshu[2]<-1
  senshu[3]<-1
  
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
  box(lwd=2.0)
  lines(c(1,6),c(-2.0,-2.0),lwd=2,lty="dashed",col="gray")
  text(6,-6,name.res[i])

  if (i == 8 ) {
    xrange.axis <- c(1,8,7)
    axis(1,xaxp=xrange.axis,lwd=2.0)
  }
  if (i== 9 || i== 11 ) {
    xrange.axis <- c(0,8,8)
    axis(1,xaxp=xrange.axis,lwd=2.0)
  }
  if (i == 10 || i == 12 ) {
    axis(1,xaxp=xrange,lwd=2.0)
  }
  
  if (i==1 || i==6 || i==11 || i==16 || i==17 ) {
    axis(2,yaxp=yrange.axis,lwd=2.0)
  }
}

name.out=paste("~/gakkai/seibutu_kanto/figs2.econdependnce_dt=1-8fs_wtune_dihed_lambda=0.1_T=300K_commp_14dipep_2012-03-01.R",sep='')
OutTiff(name.out,w_val=750,h_val=500)

