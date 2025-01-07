 source("~/Rspa/ini.R")
 source("~/Rspa/set_enecon.sh")
 source("~/Rspa/set_res.R")
 source("~/Rspa/plmGeneral_wsrange_logu.R")

 leg.pos="topright"
 is.leg <- 0

 T <- "600"
 temp <- "300"

 xrange <- c(0,14,14)
 xrange.axis <- c(0,14,14)
 yrange <- c(-14,0,14) #c(-8,0,8)
 yrange.axis <- yrange

 num.res=4
 numini=10
 numsim=2

 direconfree <- c(	"/home/yamamori/calspa/TAMD_econ/dipep_2014-07-01_ABAMD_free_end_fixed_end/ALA/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free", "/home/yamamori/calspa/TAMD_econ/dipep_2014-07-01_ABAMD_free_end_fixed_end/LEU/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free", "/home/yamamori/calspa/TAMD_econ/dipep_2014-07-01_ABAMD_free_end_fixed_end/ILE/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free", "/home/yamamori/calspa/TAMD_econ/dipep_2014-07-01_ABAMD_free_end_fixed_end/GLN/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free"  )
 direconfixed <- c(	"/home/yamamori/calspa/TAMD_econ/dipep_2014-07-01_ABAMD_free_end_fixed_end/ALA/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-07-01_ABAMD_free_end_fixed_end/LEU/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-07-01_ABAMD_free_end_fixed_end/ILE/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-07-01_ABAMD_free_end_fixed_end/GLN/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_Term_fixed"  )
 direconfixed2 <- c(	"/home/yamamori/calspa/TAMD_econ/dipep_2014-07-01_ABAMD_free_end_fixed_end/ALA/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_", "/home/yamamori/calspa/TAMD_econ/dipep_2014-07-01_ABAMD_free_end_fixed_end/LEU/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_", "/home/yamamori/calspa/TAMD_econ/dipep_2014-07-01_ABAMD_free_end_fixed_end/ILE/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_", "/home/yamamori/calspa/TAMD_econ/dipep_2014-07-01_ABAMD_free_end_fixed_end/GLN/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_NH2_Term_fixed"  )
 proname <- c(	"ALADv", "LEUDv", "ILEDv", "GLNDv"  )

 file.name=paste("/home/yamamori/seminars/GS/2014-07-01//eps/fig_econdependence_ABAMD_4AA_Term_fixed_vs_woTerm_fixed_2014-07-01.eps",sep='')
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

  file.names[1]=paste(direconfree[i],"/",proname[i],"_econ_300_100ps_1-",numini,".econ.rmsd.av",sep='')
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
  text(8,-12,proname[i],cex=1.0)
  box(lwd=2.0)

  lines(c(1,6),c(-1.0,-1.0),lwd=2,lty="dashed",col="gray")
  lines(c(7,7),c(-10,0),lwd=2,lty="dashed",col="gray")
  lines(c(8,8),c(-10,0),lwd=2,lty="dashed",col="gray")
#  lines(c(10,10),c(-10,0),lwd=2,lty="dashed",col="gray")

  if ( i > -1 ) {
    axis(1,xaxp=xrange.axis,lwd=2.0)
  }

  if ( i %% 4 == 1 || i == 1) {
    axis(2,yaxp=yrange.axis,lwd=2.0)
  }
    
}

mtext(outer=T,label.x,side=1,line=2.0,cex=1.0)
mtext(outer=T,label.y,side=2,line=2.0,cex=1.0)

