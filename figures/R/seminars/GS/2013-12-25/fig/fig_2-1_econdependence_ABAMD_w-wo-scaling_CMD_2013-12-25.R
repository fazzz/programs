source("~/Rspa/ini.R")

source("~/Rspa/set_enecon.sh")
source("~/Rspa/set_res.R")
source("~/Rspa/plmGeneral_wsrange.R")
source("~/Rspa/plmGeneral_wsrange_log.R")

leg.pos="topright"
is.leg <- 0

temp <- "300"

xrange <- c(0,13,12)
yrange <- c(-20,1,21)

name.res <- c( "ASN", "GLN", "VAL", "ILE" )

num.res <- length(name.res)

numini=10

dirbaseABAwscaling  <- "~/calspa/TAMD_econ/dipep_2013-12-01_sander_af_minimize_all/"
dirbaseABAwoscaling <- "~/calspa/TAMD_econ/dipep_2013-12-01_sander_af_minimize_all/"
dirbaseCMD          <- "~/calspa/refcalc/dipep/"

file.name=paste("~/seminars/GS/2013-12-25/eps/fig_2-1_econdependence_ABAMD_w-wo-scaling_CMD_2013-12-25.eps",sep='')
#postscript(file.name,width=8,height=6,horizontal=FALSE,onefile=FALSE,paper="special")
postscript(file.name,width=8,height=3,horizontal=FALSE,onefile=FALSE,paper="special")

hutosa <- 2

sen <- NULL
iro <- NULL
id.ys <- NULL
ids.ys <- NULL
tenshu <- NULL
senshu <- NULL

num.line <- 3

for ( k in 1:num.line  ) {
  sen[k]<-c(1)
  id.ys[k]<-c(2)
  ids.ys[k]<-c(3)
  tenshu[k]<-c(20)
}

iro[1]<-"red"
iro[2]<-"red"
iro[3]<-"black"

senshu[1]<-1
senshu[2]<-2
senshu[3]<-1

tune <- c(0.1)
ntune <- length(tune)

ntc<-"2"
T<-"600"

ffname=paste("Amber_tuned",tune,sep='')

par(mfrow=c(1,4))
#par(oma=c(7.5,4.0,2.0,2.0))
par(oma=c(2.0,4.0,2.0,2.0))

file.names<-NULL
for (i in 1:num.res) {

  par(mar=c(0.0,0.0,0.0,0.0))

#  file.names[1]=paste(dirbaseABAwscaling,name.res[i],"/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wff",ffname,"_2013-12-01/",name.res[i],"Dv_econ_",temp[1],"_100ps_sa=",T,"_1-",numini,".econ_2.rmsd.av",sep='')
#  file.names[2]=paste(dirbaseABAwoscaling,name.res[i],"/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_2013-12-02/",name.res[i],"Dv_econ_",temp[1],"_100ps_sa=",T,"_1-",numini,".econ_2.rmsd.av",sep='')
#  file.names[3]=paste(dirbaseCMD,name.res[i],"/econ_wdi_NVE_100ps_ECsander_SAsander_afclust_and_minimize_all_2013-12-12/",name.res[i],"Dv_econ_",temp[1],"_100ps_1-",numini,"_",ntc,"ntc.econ_2.rmsd.av",sep='')

  file.names[1]=paste(dirbaseABAwscaling,name.res[i],"/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wff",ffname,"_2013-12-01/",name.res[i],"Dv_econ_",temp[1],"_100ps_sa=",T,"_1-",numini,".econ.rmsd.av",sep='')
  file.names[2]=paste(dirbaseABAwoscaling,name.res[i],"/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_2013-12-02/",name.res[i],"Dv_econ_",temp[1],"_100ps_sa=",T,"_1-",numini,".econ.rmsd.av",sep='')
  file.names[3]=paste(dirbaseCMD,name.res[i],"/econ_wdi_NVE_100ps_ECsander_SAsander_afclust_and_minimize_all_2013-12-12/",name.res[i],"Dv_econ_",temp[1],"_100ps_1-",numini,"_",ntc,"ntc.econ.rmsd.av",sep='')

  plmGeneralwsrange(data.names=file.names,
                    sd.names=file.names,
                    id.ys=id.ys,
                    ids.ys=ids.ys,
                    label.size=0.5,axis.size=2.0,
                    iro=iro,axis.ft="F",is.header="T",
                    sdiro=iro,
                    xrange=xrange,
                    yrange=yrange,
                    sdyrange=yrange,
                    is.sen=sen,width=10.0,
                    warrow="T")
  text(6,-6.0,name.res[i],cex=1.5)
  box(lwd=2.0)

  lines(c(1,6),c(-1.0,-1.0),lwd=2,lty="dashed",col="gray")

  xrange.axis <- c(0,12,12)
  axis(1,xaxp=xrange.axis,lwd=2.0)

  if (i==1 ) {
    yrange.axis <- c(-20,0,20)
    axis(2,yaxp=yrange.axis,lwd=2.0)
  }
    
}
