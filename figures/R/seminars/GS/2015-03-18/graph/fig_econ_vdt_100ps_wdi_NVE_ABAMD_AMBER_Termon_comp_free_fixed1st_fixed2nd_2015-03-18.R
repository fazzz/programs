 source("~/Rspa/ini.R")
 source("~/Rspa/set_enecon.sh")
 source("~/Rspa/set_res.R")
 source("~/Rspa/plmGeneral_wsrange_logu.R")

 leg.pos="topright"
 is.leg <- 0

 T <- "600"
 temp <- "300"

 xrange <- c(0,16,16)
 xrange.axis <- c(0,14,14)
 yrange <-  c(-8,0,8) #c(-14,0,14)
 yrange.axis <- c(-8,-1,7)

 num.res=27
 numini=10
 numsim=15

 direconfree <- c(	"/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/GLY/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/ALA/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/VAL/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/LEU/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/ILE/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/PRO/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/MET/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_S_free", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/CYS/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/CYX/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/CYM/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/SER/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/THR/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/ASN/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/ASH/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/GLN/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/GLH/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/ASP/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/GLU/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/LYS/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/LYN/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/ARG/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_C_free", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/HIE/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/HIP/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/HID/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/PHE/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/TY3/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/TRP/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free"  )
 direconfixed <- c(	"/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/GLY/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/ALA/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/VAL/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/LEU/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/ILE/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/PRO/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/MET/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/CYS/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/CYX/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_OH_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/CYM/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/SER/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/THR/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/ASN/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/ASH/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/GLN/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/GLH/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/ASP/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/GLU/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/LYS/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_NH3_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/LYN/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/ARG/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_C_free_CH3_fixd", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/HIE/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/HIP/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/HID/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/PHE/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/TY3/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/TRP/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_Term_fixed"  )
 direconfixed2 <- c(	"XX", "XX", "XX", "XX", "XX", "XX", "XX", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/CYS/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_SH_Term_fixed", "XX", "XX", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/SER/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_OH_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/THR/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_OH_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/ASN/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_NH2_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/ASH/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_OH_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/GLN/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_NH2_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/GLH/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_OH_Term_fixed", "XX", "XX", "XX", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/LYN/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_NH2_Term_fixed", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/ARG/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_NH2_Term_fixed", "XX", "XX", "XX", "XX", "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/TY3/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_CH3_OH_Term_fixed", "XX"  )
 proname <- c(	"GLYDv", "ALADv", "VALDv", "LEUDv", "ILEDv", "PRODv", "METDv", "CYSDv", "CYXDv", "CYMDv", "SERDv", "THRDv", "ASNDv", "ASHDv", "GLNDv", "GLHDv", "ASPDv", "GLUDv", "LYSDv", "LYNDv", "ARGDv", "HIEDv", "HIPDv", "HIDDv", "PHEDv", "TY3Dv", "TRPDv"  )
 resname <- c(	"GLY", "ALA", "VAL", "LEU", "ILE", "PRO", "MET", "CYS", "CYX", "CYM", "SER", "THR", "ASN", "ASH", "GLN", "GLH", "ASP", "GLU", "LYS", "LYN", "ARG", "HIE", "HIP", "HID", "PHE", "TY3", "TRP"  )

 file.name=paste("/home/yamamori/seminars/GS/2015-03-18/eps/fig_econ_vdt_100ps_wdi_NVE_ABAMD_AMBER_Termon_comp_free_fixed1st_fixed2nd_2015-03-18.eps",sep='')
 postscript(file.name,width=14.45,height=5.50,horizontal=FALSE,onefile=FALSE,paper="special")

 label.x <- "dt (fs)"
 label.y <- "log(dE)"

 par(mfrow=c(3,9))
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
  text(8,-7,resname[i],cex=1.75)
  box(lwd=2.0)

#  lines(c(1,6),c(-1.0,-1.0),lwd=2,lty="dashed",col="gray")
#  lines(c(7,7),c(-10,0),lwd=2,lty="dashed",col="gray")
#  lines(c(8,8),c(-10,0),lwd=2,lty="dashed",col="gray")
#  lines(c(10,10),c(-10,0),lwd=2,lty="dashed",col="gray")

  if ( i == 19 ) {
    xrange.axis <- c(0,15,15)
    axis(1,xaxp=xrange.axis,lwd=2.0)
  }
  else if ( i > 18 && i < 27 ) {
    xrange.axis <- c(1,16,15)
    axis(1,xaxp=xrange.axis,lwd=2.0)
  }
  else if ( i > 18 ) {
    xrange.axis <- c(1,16,15)
    axis(1,xaxp=xrange.axis,lwd=2.0)
  }

  if ( i == 1) {
    yrange.axis <- c(-8,2,10)
    axis(2,yaxp=yrange.axis,lwd=2.0)
  }
  else if ( i == 10) {
    yrange.axis <- c(-8,1,9)
    axis(2,yaxp=yrange.axis,lwd=2.0)
  }
  else if ( i == 19) {
    yrange.axis <- c(-8,1,9)
    axis(2,yaxp=yrange.axis,lwd=2.0)
  }
    
}

#mtext(outer=T,label.x,side=1,line=4.0,cex=1.0)
#mtext(outer=T,label.y,side=2,line=3.0,cex=1.0)

