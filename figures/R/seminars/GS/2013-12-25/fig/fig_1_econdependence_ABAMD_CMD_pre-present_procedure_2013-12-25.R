source("~/Rspa/ini.R")

source("~/Rspa/set_enecon.sh")
source("~/Rspa/set_res.R")
source("~/Rspa/plmGeneral_wsrange.R")
source("~/Rspa/plmGeneral_wsrange_log.R")

leg.pos="topright"
is.leg <- 0

temp <- "300"

xrange <- c(0,13,13)
#yrange <- c(-7,0,7)
yrange <- c(-20,0,20)

name.res <- c( "ALA", "ASN", "GLY", "GLN")

#name.res <- c( "ALA" )

num.res <- length(name.res)

numini=10

dirbaseABApept <- "~/calspa/TAMD_econ/dipep_2013-12-01_sander_af_minimize_all/"
dirbaseABApept2 <- "/misc/ivy3/yamamori/ajisai/TAMD_econ/TAMD_econ/dipep/"

file.name=paste("~/seminars/GS/2013-12-25/eps/fig_1_econdependence_ABAMD_CMD_pre-present_procedure_2013-12-25.eps",sep='')
postscript(file.name,width=5.6,height=4.2,horizontal=FALSE,onefile=FALSE,paper="special")

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
iro[2]<-"black"

senshu[1]<-1
senshu[2]<-1

tune <- c(0.1)
ntune <- length(tune)

T<-"600"

ffname=paste("Amber_tuned",tune,sep='')

#par(mfrow=c(4,5))
par(mfrow=c(2,2))
par(oma=c(2.0,4.0,2.0,2.0))

file.names<-NULL
for (i in 1:num.res) {

  par(mar=c(0.0,0.0,0.0,0.0))

#  file.names[1]=paste(dirbaseABApept,name.res[i],"/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wff",ffname,"_2013-12-01/",name.res[i],"Dv_econ_",temp[1],"_100ps_sa=",T,"_1-",numini,".econ_2.rmsd.av",sep='')
#  file.names[2]=paste(dirbaseABApept2,name.res[i],"/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_2012-04-17_pc6/dtune=",tune[1],"/",name.res[i],"Dv_econ_",temp[1],"_100ps_1-",numini,".econ_2.rmsd.av",sep='')

  file.names[1]=paste(dirbaseABApept,name.res[i],"/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wff",ffname,"_2013-12-01/",name.res[i],"Dv_econ_",temp[1],"_100ps_sa=",T,"_1-",numini,".econ.rmsd.av",sep='')
  file.names[2]=paste(dirbaseABApept2,name.res[i],"/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_2012-04-17_pc6/dtune=",tune[1],"/",name.res[i],"Dv_econ_",temp[1],"_100ps_1-",numini,".econ.rmsd.av",sep='')

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

  lines(c(1,14),c(-1.0,-1.0),lwd=2,lty="dashed",col="gray")

#  if (i==16 || i==17 || i==18 || i==19 ) {
  if (i==3 ) {
    xrange.axis <- c(0,13,13)
    axis(1,xaxp=xrange.axis,lwd=2.0)
  }

#  if (i==20) {
  if (i==4) {
    xrange.axis <- c(0,13,13)
    axis(1,xaxp=xrange.axis,lwd=2.0)
  }

  if (i==1 ) {
#    yrange.axis <- c(-7,0,7)
    yrange.axis <- c(-18,0,18)
    axis(2,yaxp=yrange.axis,lwd=2.0)
  }
  
#  if (i==6 || i==11 || i==16 ) {
  if (i==3  ) {
#    yrange.axis <- c(-7,-1,6)
    yrange.axis <- c(-20,0,20)
    axis(2,yaxp=yrange.axis,lwd=2.0)
  }
    
}
