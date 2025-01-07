 source("~/Rspa/ini.R")
 source("~/Rspa/set_enecon.sh")
 source("~/Rspa/set_res.R")
 source("~/Rspa/plmGeneral_wsrange_logu.R")

 leg.pos="topright"
 is.leg <- 0

 T <- "600"
 temp <- "300"

 xrange <- c(0,14,14)
 xrange.axis <- c(0,13,13)
 xrange.axist <- c(1,13,12)
 yrange <- c(-9,0,9)
 yrange.axis <- c(-8,0,8)

 num.res=4
 numini=10
 numsim=16

 direcon <- c(	"/home/yamamori/calspa/TAMD_econ/dipep_2014-03-05_sander_af_minimize/GLY/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning", "/home/yamamori/calspa/TAMD_econ/dipep_2014-03-05_sander_af_minimize/ALA/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning", "/home/yamamori/calspa/TAMD_econ/dipep_2014-03-05_sander_af_minimize/ASP/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning", "/home/yamamori/calspa/TAMD_econ/dipep_2014-03-05_sander_af_minimize/GLU/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning"  )
 direconfixed <- c(	"/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/GLY/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/ALA/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/ASP/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/GLU/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_fixed"  )
 direconfixed2 <- c(	"/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/GLY/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_XX", "/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/ALA/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_XX", "/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/ASP/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_XX", "/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/GLU/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_XX"  )
 proname <- c(	"GLYDv", "ALADv", "ASPDv", "GLUDv"  )
 pronamet <- c(	"GLY", "ALA", "ASP", "GLU"  )

 file.name=paste("/home/yamamori/seminars/GS/2014-09-17/eps/fig_econservation_dependence_ABAMD_free_terminal_vs_1st_2nd_fixed_terminal_before_debig_2014-09-17.eps",sep='')
 postscript(file.name,width=6.95,height=2.50,horizontal=FALSE,onefile=FALSE,paper="special")

 label.x <- "dt (fs)"
 label.y <- "log(dE)"

 par(mfrow=c(1,4))
 par(omi=c(0.75,0.80,0.2,0.2))

 hutosa <- c(2, 2, 2)
 sen <- c(3, 3, 3)
 id.ys  <- c(2, 2, 2)
 ids.ys <- c(3, 3, 3)
 tenshu <- c(20, 20, 20)
 senshu <- c(1, 1, 1)
 iro    <-c( "red", "black", "green" )
 file.names<-NULL

 for (i in 1:num.res) {
  par(mar=c(0.0,0.0,0.0,0.0))

  file.names[1]=paste(direcon[i],"/",proname[i],"_econ_300_100ps_1-",numini,".econ.rmsd.av",sep='')
  file.names[2]=paste(direconfixed[i],"/",proname[i],"_econ_300_100ps_1-",numini,".econ.rmsd.av",sep='')
  file.names[3]=paste(direconfixed2[i],"/",proname[i],"_econ_300_100ps_1-",numini,".econ.rmsd.av",sep='')

  plmGeneralwsrangelogu(data.names=file.names,
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
  text(8,-6,pronamet[i],cex=1.75)
  box(lwd=2.0)

  lines(c(1,6),c(-1.0,-1.0),lwd=2,lty="dashed",col="gray")
  lines(c(7,7),c(-10,0),lwd=2,lty="dashed",col="gray")
  lines(c(8,8),c(-10,0),lwd=2,lty="dashed",col="gray")
#  lines(c(10,10),c(-10,0),lwd=2,lty="dashed",col="gray")

  if ( i == -1 ) {
    axis(1,xaxp=xrange.axist,lwd=2.0)
  }

  if ( i > -1 ) {
    axis(1,xaxp=xrange.axis,lwd=2.0)
  }

  if ( i %% 4 == 1 || i == 1) {
    axis(2,yaxp=yrange.axis,lwd=2.0)
  }
    
}

mtext(outer=T,label.x,side=1,line=4.0,cex=1.0)
mtext(outer=T,label.y,side=2,line=3.0,cex=1.0)

