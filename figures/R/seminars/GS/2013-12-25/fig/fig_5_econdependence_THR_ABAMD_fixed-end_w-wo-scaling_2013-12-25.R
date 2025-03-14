source("~/Rspa/ini.R")

source("~/Rspa/set_enecon.sh")
source("~/Rspa/set_res.R")
source("~/Rspa/plmGeneral_wsrange.R")
source("~/Rspa/plmGeneral_wsrange_log.R")

leg.pos="topright"
is.leg <- 0

T <- "600"
temp <- "300"

clustname.CH3 <- "_CH3"
clustname.CH3.OH <- "_CH3_OH"

xrange <- c(0,13,12)
yrange <- c(-20,1,21)

name.res <- c( "THR" )

num.res <- length(name.res)

numini=10

dirbaseABApept.wscaling    <- "~/calspa/TAMD_econ/dipep_2013-12-01_sander_af_minimize_all/"
dirbaseABApept.fixedend  <- "~/calspa/TAMD_econ/dipep_2013-12-16_sander_af_minimize_all_fixed_end/"
dirbaseABApept.woscaling <- "~/calspa/TAMD_econ/dipep_2013-12-01_sander_af_minimize_all/"

file.name=paste("~/seminars/GS/2013-12-25/eps/fig_5_econdependence_THR_ABAMD_fixed-end_w-wo-scaling_2013-12-25.eps",sep='')
postscript(file.name,width=5,height=4,horizontal=FALSE,onefile=FALSE,paper="special")

hutosa <- 2

num.line <- 4

sen <- NULL
id.ys <- NULL
ids.ys <- NULL
tenshu <- NULL

for ( k in 1:num.line  ) {
  sen[k]<-c(1)
  id.ys[k]<-c(2)
  ids.ys[k]<-c(3)
  tenshu[k]<-c(20)
}

iro[1]<-"blue"
iro[2]<-"blue"
iro[3]<-"red"
iro[4]<-"red"

senshu[1]<-1
senshu[2]<-2
senshu[3]<-1
senshu[4]<-2

tune <- c(0.1)
ntune <- length(tune)

T<-"600"

ffname=paste("Amber_tuned",tune,sep='')

par(mfrow=c(1,1))
par(oma=c(2.0,4.0,2.0,2.0))

file.names<-NULL
for (i in 1:num.res) {

  par(mar=c(0.0,0.0,0.0,0.0))

#  file.names[1]=paste(dirbaseABApept.fixedend,name.res[i],"/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_2013-12-16_",clustname.CH3,"/",name.res[i],"Dv_econ_",temp[1],"_100ps_sa=",T,"_1-",numini,".econ_2.rmsd.av",sep='')
#  file.names[2]=paste(dirbaseABApept.fixedend,name.res[i],"/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_2013-12-16_",clustname.CH3.OH,"/",name.res[i],"Dv_econ_",temp[1],"_100ps_sa=",T,"_1-",numini,".econ_2.rmsd.av",sep='')
#  file.names[3]=paste(dirbaseABApept.wscaling,name.res[i],"/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wff",ffname,"_2013-12-01/",name.res[i],"Dv_econ_",temp[1],"_100ps_sa=",T,"_1-",numini,".econ_2.rmsd.av",sep='')
#  file.names[4]=paste(dirbaseABApept.woscaling,name.res[i],"/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_2013-12-02/",name.res[i],"Dv_econ_",temp[1],"_100ps_sa=",T,"_1-",numini,".econ_2.rmsd.av",sep='')

  file.names[1]=paste(dirbaseABApept.fixedend,name.res[i],"/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_2013-12-16_",clustname.CH3,"/",name.res[i],"Dv_econ_",temp[1],"_100ps_sa=",T,"_1-",numini,".econ.rmsd.av",sep='')
  file.names[2]=paste(dirbaseABApept.fixedend,name.res[i],"/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_2013-12-16_",clustname.CH3.OH,"/",name.res[i],"Dv_econ_",temp[1],"_100ps_sa=",T,"_1-",numini,".econ.rmsd.av",sep='')
  file.names[3]=paste(dirbaseABApept.wscaling,name.res[i],"/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wff",ffname,"_2013-12-01/",name.res[i],"Dv_econ_",temp[1],"_100ps_sa=",T,"_1-",numini,".econ.rmsd.av",sep='')
  file.names[4]=paste(dirbaseABApept.woscaling,name.res[i],"/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_2013-12-02/",name.res[i],"Dv_econ_",temp[1],"_100ps_sa=",T,"_1-",numini,".econ.rmsd.av",sep='')

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
  text(10,-16.0,name.res[i],cex=1.5)
  box(lwd=2.0)

  lines(c(1,6),c(-1.0,-1.0),lwd=2,lty="dashed",col="gray")

  xrange.axis <- c(0,12,12)
  axis(1,xaxp=xrange.axis,lwd=2.0)

  yrange.axis <- c(-20,0,20)
  axis(2,yaxp=yrange.axis,lwd=2.0)

}
