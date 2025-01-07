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
 yrange <- c(-8,0,8)
 yrange.axis <- yrange

 num.res=20
 numini=10
 numsim=16

 direcon <- c(	"/home/yamamori/calspa/TAMD_econ/dipep_2014-03-05_sander_af_minimize/ALA/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning", "/home/yamamori/calspa/TAMD_econ/dipep_2014-03-05_sander_af_minimize/VAL/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning", "/home/yamamori/calspa/TAMD_econ/dipep_2014-03-05_sander_af_minimize/LEU/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning", "/home/yamamori/calspa/TAMD_econ/dipep_2014-03-05_sander_af_minimize/ILE/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning", "/home/yamamori/calspa/TAMD_econ/dipep_2014-03-05_sander_af_minimize/THR/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning", "/home/yamamori/calspa/TAMD_econ/dipep_2014-03-05_sander_af_minimize/SER/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning", "/home/yamamori/calspa/TAMD_econ/dipep_2014-03-05_sander_af_minimize/ASN/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning", "/home/yamamori/calspa/TAMD_econ/dipep_2014-03-05_sander_af_minimize/GLN/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning", "/home/yamamori/calspa/TAMD_econ/dipep_2014-03-05_sander_af_minimize/LYS/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning", "/home/yamamori/calspa/TAMD_econ/dipep_2014-03-05_sander_af_minimize/CYS/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning", "/home/yamamori/calspa/TAMD_econ/dipep_2014-03-05_sander_af_minimize/ASP/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning", "/home/yamamori/calspa/TAMD_econ/dipep_2014-03-05_sander_af_minimize/ARG/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning", "/home/yamamori/calspa/TAMD_econ/dipep_2014-03-05_sander_af_minimize/GLU/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning", "/home/yamamori/calspa/TAMD_econ/dipep_2014-03-05_sander_af_minimize/GLY/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning", "/home/yamamori/calspa/TAMD_econ/dipep_2014-03-05_sander_af_minimize/HIE/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning", "/home/yamamori/calspa/TAMD_econ/dipep_2014-03-05_sander_af_minimize/MET/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning", "/home/yamamori/calspa/TAMD_econ/dipep_2014-03-05_sander_af_minimize/TRP/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning", "/home/yamamori/calspa/TAMD_econ/dipep_2014-03-05_sander_af_minimize/TYR/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_free_OH", "/home/yamamori/calspa/TAMD_econ/dipep_2014-03-05_sander_af_minimize/PHE/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning", "/home/yamamori/calspa/TAMD_econ/dipep_2014-03-05_sander_af_minimize/PRO/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_2014-03-17"  )
 direconfixed <- c(	"/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/ALA/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/VAL/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/LEU/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/ILE/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/THR/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/SER/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/ASN/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_NH2_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/GLN/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_NH2_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/LYS/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_NH3_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/CYS/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_SH_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/ASP/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/ARG/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/GLU/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/GLY/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/HIE/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/MET/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/TRP/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/TYR/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_XX", "/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/PHE/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/PRO/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_fixed_2014-03-17"  )
 direconfixed2 <- c(	"/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/ALA/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_XX", "/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/VAL/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_XX", "/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/LEU/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_XX", "/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/ILE/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_XX", "/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/THR/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_OH_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/SER/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_OH_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/ASN/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_XX", "/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/GLN/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_XX", "/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/LYS/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_XX", "/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/CYS/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_XX", "/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/ASP/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_XX", "/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/ARG/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_XX", "/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/GLU/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_XX", "/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/GLY/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_XX", "/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/HIE/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_XX", "/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/MET/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_XX", "/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/TRP/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_XX", "/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/TYR/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/PHE/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_XX", "/home/yamamori/calspa/TAMD_econ/dipep_2014-02-25_sander_af_minimize_all_fixed_end/PRO/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_2014-03-17"  )
 proname <- c(	"ALADv", "VALDv", "LEUDv", "ILEDv", "THRDv", "SERDv", "ASNDv", "GLNDv", "LYSDv", "CYSDv", "ASPDv", "ARGDv", "GLUDv", "GLYDv", "HIEDv", "METDv", "TRPDv", "TYRDv", "PHEDv", "PRODv"  )

 file.name=paste("/home/yamamori/seminars/GS/2014-04-08//eps/fig_econdependence_fixed_Term_vs_free_Term_20-dipeptides_2014-04-08.eps",sep='')
 postscript(file.name,width=8.45,height=7.00,horizontal=FALSE,onefile=FALSE,paper="special")

 label.x <- "dt (fs)"
 label.y <- "log(dE)"

 par(mfrow=c(4,5))
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
  text(8,-6,proname[i],cex=1.0)
  box(lwd=2.0)

  lines(c(1,6),c(-1.0,-1.0),lwd=2,lty="dashed",col="gray")
  lines(c(7,7),c(-10,0),lwd=2,lty="dashed",col="gray")
  lines(c(8,8),c(-10,0),lwd=2,lty="dashed",col="gray")
#  lines(c(10,10),c(-10,0),lwd=2,lty="dashed",col="gray")

  if ( i > 15 ) {
    axis(1,xaxp=xrange.axis,lwd=2.0)
  }

  yrange.axis <- c(-7,0,7)
  if ( i == 1 ) {
    axis(2,yaxp=yrange.axis,lwd=2.0)
  }
  yrange.axis <- c(-8,-1,7)
  if (i == 6 || i == 11) {
    axis(2,yaxp=yrange.axis,lwd=2.0)
  }
  yrange.axis <- c(-8,-1,7)
  if ( i == 16) {
    axis(2,yaxp=yrange.axis,lwd=2.0)
  }
    
}

mtext(outer=T,label.x,side=1,line=2.0,cex=1.0)
mtext(outer=T,label.y,side=2,line=2.0,cex=1.0)

