 source("~/Rspa/ini.R")
 source("~/Rspa/set_enecon.sh")
 source("~/Rspa/set_res.R")
 source("~/Rspa/plmGeneral_wsrange_logu.R")

 leg.pos="topright"
 is.leg <- 0

 T <- "600"
 temp <- "300"

 xrange <- c(0,16,16)
 xrange.axis <- c(0,16,16)
# yrange <-  c(-8,0,8) #c(-14,0,14)
 yrange <-  c(-8,2,10) #c(-14,0,14)
 yrange.axis <- yrange

 num.res=27
 numini=10
 numsim=15

 direcon <- c(	"/home/yamamori/calspa/TAMD_econ/dipep_2015-01-28_ABAMD_dubugged_mod_parm//GLY/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_PC6_tune_s_0.1", "/home/yamamori/calspa/TAMD_econ/dipep_2015-01-28_ABAMD_dubugged_mod_parm//ALA/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_PC6_tune_s_0.1", "/home/yamamori/calspa/TAMD_econ/dipep_2015-01-28_ABAMD_dubugged_mod_parm//VAL/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_PC6_tune_s_0.1", "/home/yamamori/calspa/TAMD_econ/dipep_2015-01-28_ABAMD_dubugged_mod_parm//LEU/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_PC6_tune_s_0.1", "/home/yamamori/calspa/TAMD_econ/dipep_2015-01-28_ABAMD_dubugged_mod_parm//ILE/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_PC6_tune_s_0.1", "/home/yamamori/calspa/TAMD_econ/dipep_2015-01-28_ABAMD_dubugged_mod_parm//PRO/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_PC6_tune_s_0.1", "/home/yamamori/calspa/TAMD_econ/dipep_2015-01-28_ABAMD_dubugged_mod_parm//MET/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_PC6_tune_s_0.1", "/home/yamamori/calspa/TAMD_econ/dipep_2015-01-28_ABAMD_dubugged_mod_parm//CYS/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_PC6_tune_s_0.1", "/home/yamamori/calspa/TAMD_econ/dipep_2015-01-28_ABAMD_dubugged_mod_parm//CYX/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_PC6_tune_s_0.1", "/home/yamamori/calspa/TAMD_econ/dipep_2015-01-28_ABAMD_dubugged_mod_parm//CYM/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_PC6_tune_s_0.1", "/home/yamamori/calspa/TAMD_econ/dipep_2015-01-28_ABAMD_dubugged_mod_parm//SER/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_PC6_tune_s_0.1", "/home/yamamori/calspa/TAMD_econ/dipep_2015-01-28_ABAMD_dubugged_mod_parm//THR/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_PC6_tune_s_0.1", "/home/yamamori/calspa/TAMD_econ/dipep_2015-01-28_ABAMD_dubugged_mod_parm//ASN/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_PC6_tune_s_0.1", "/home/yamamori/calspa/TAMD_econ/dipep_2015-01-28_ABAMD_dubugged_mod_parm//ASH/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_PC6_tune_s_0.1", "/home/yamamori/calspa/TAMD_econ/dipep_2015-01-28_ABAMD_dubugged_mod_parm//GLN/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_PC6_tune_s_0.1", "/home/yamamori/calspa/TAMD_econ/dipep_2015-01-28_ABAMD_dubugged_mod_parm//GLH/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_PC6_tune_s_0.1", "/home/yamamori/calspa/TAMD_econ/dipep_2015-01-28_ABAMD_dubugged_mod_parm//ASP/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_PC6_tune_s_0.1", "/home/yamamori/calspa/TAMD_econ/dipep_2015-01-28_ABAMD_dubugged_mod_parm//GLU/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_PC6_tune_s_0.1", "/home/yamamori/calspa/TAMD_econ/dipep_2015-01-28_ABAMD_dubugged_mod_parm//LYS/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_PC6_tune_s_0.1", "/home/yamamori/calspa/TAMD_econ/dipep_2015-01-28_ABAMD_dubugged_mod_parm//LYN/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_PC6_tune_s_0.1", "/home/yamamori/calspa/TAMD_econ/dipep_2015-01-28_ABAMD_dubugged_mod_parm//ARG/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_PC6_tune_s_0.1", "/home/yamamori/calspa/TAMD_econ/dipep_2015-01-28_ABAMD_dubugged_mod_parm//HIE/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_PC6_tune_s_0.1", "/home/yamamori/calspa/TAMD_econ/dipep_2015-01-28_ABAMD_dubugged_mod_parm//HIP/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_PC6_tune_s_0.1", "/home/yamamori/calspa/TAMD_econ/dipep_2015-01-28_ABAMD_dubugged_mod_parm//HID/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_PC6_tune_s_0.1", "/home/yamamori/calspa/TAMD_econ/dipep_2015-01-28_ABAMD_dubugged_mod_parm//PHE/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_PC6_tune_s_0.1", "/home/yamamori/calspa/TAMD_econ/dipep_2015-01-28_ABAMD_dubugged_mod_parm//TY3/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_PC6_tune_s_0.1", "/home/yamamori/calspa/TAMD_econ/dipep_2015-01-28_ABAMD_dubugged_mod_parm//TRP/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_PC6_tune_s_0.1"  )
 proname <- c(	"GLYDv", "ALADv", "VALDv", "LEUDv", "ILEDv", "PRODv", "METDv", "CYSDv", "CYXDv", "CYMDv", "SERDv", "THRDv", "ASNDv", "ASHDv", "GLNDv", "GLHDv", "ASPDv", "GLUDv", "LYSDv", "LYNDv", "ARGDv", "HIEDv", "HIPDv", "HIDDv", "PHEDv", "TY3Dv", "TRPDv"  )
 resname <- c(	"GLY", "ALA", "VAL", "LEU", "ILE", "PRO", "MET", "CYS", "CYX", "CYM", "SER", "THR", "ASN", "ASH", "GLN", "GLH", "ASP", "GLU", "LYS", "LYN", "ARG", "HIE", "HIP", "HID", "PHE", "TY3", "TRP"  )

 file.name=paste("/home/yamamori/seminars/GS/2015-03-18/eps/fig_econ_vdt_100ps_wdi_NVE_ABAMD_AMBER_Termon_comp_woscale_wscale_2015-03-18.eps",sep='')
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
 iro    <-c( "red", "blue", "green" )
 file.names<-NULL

 for (i in 1:num.res) {
  par(mar=c(0.0,0.0,0.0,0.0))

  file.names[1]=paste(direcon[i],"/",proname[i],"_econ_300_100ps_1-",numini,".econ.rmsd.av",sep='')

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

