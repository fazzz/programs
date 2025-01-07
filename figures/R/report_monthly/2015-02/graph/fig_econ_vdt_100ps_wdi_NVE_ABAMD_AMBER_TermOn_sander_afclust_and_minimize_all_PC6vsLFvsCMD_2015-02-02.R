 source("~/Rspa/ini.R")
 source("~/Rspa/set_enecon.sh")
 source("~/Rspa/set_res.R")
 source("~/Rspa/plmGeneral_wsrange_logu.R")

 leg.pos="topright"
 is.leg <- 0

 xrange <- c(0,16,16)
 yrange <- c(-8,2,10)

 num.res=27

 RES <- c( "GLY" , "ALA" , "VAL" , "LEU" , "ILE" , "PRO" , "MET" , "CYS" , "CYX" , "CYM" , "SER" , "THR" , "ASN" , "ASH" , "GLN" , "GLH" , "ASP" , "GLU" , "LYS" , "LYN" , "ARG" , "HIE" , "HIP" , "HID" , "PHE" , "TY3" , "TRP"  )
 econPC6 <- c(	"/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/GLY/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free/GLYDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/ALA/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free/ALADv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/VAL/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_fixed/VALDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/LEU/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free/LEUDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/ILE/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free/ILEDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/PRO/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free/PRODv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/MET/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_S_free/METDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/CYS/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free/CYSDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/CYX/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free/CYXDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/CYM/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free/CYMDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/SER/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free/SERDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/THR/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free/THRDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/ASN/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free/ASNDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/ASH/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free/ASHDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/GLN/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free/GLNDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/GLH/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free/GLHDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/ASP/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free/ASPDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/GLU/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free/GLUDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/LYS/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free/LYSDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/LYN/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free/LYNDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/ARG/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_C_free/ARGDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/HIE/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free/HIEDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/HIP/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free/HIPDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/HID/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free/HIDDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/PHE/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free/PHEDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/TY3/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free/TY3Dv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/TRP/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_Term_free/TRPDv_econ_300_100ps_1-10.econ.rmsd.av"  )
 econmLF <- c(	"/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/GLY/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_LeapFrog_Term_free/GLYDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/ALA/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_LeapFrog_Term_free/ALADv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/VAL/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_LeapFrog_Term_fixed/VALDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/LEU/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_LeapFrog_Term_free/LEUDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/ILE/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_LeapFrog_Term_free/ILEDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/PRO/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_LeapFrog_Term_free/PRODv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/MET/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_LeapFrog_Term_S_free/METDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/CYS/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_LeapFrog_Term_free/CYSDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/CYX/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_LeapFrog_Term_free/CYXDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/CYM/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_LeapFrog_Term_free/CYMDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/SER/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_LeapFrog_Term_free/SERDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/THR/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_LeapFrog_Term_free/THRDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/ASN/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_LeapFrog_Term_free/ASNDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/ASH/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_LeapFrog_Term_free/ASHDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/GLN/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_LeapFrog_Term_free/GLNDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/GLH/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_LeapFrog_Term_free/GLHDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/ASP/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_LeapFrog_Term_free/ASPDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/GLU/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_LeapFrog_Term_free/GLUDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/LYS/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_LeapFrog_Term_free/LYSDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/LYN/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_LeapFrog_Term_free/LYNDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/ARG/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_LeapFrog_Term_C_free/ARGDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/HIE/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_LeapFrog_Term_free/HIEDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/HIP/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_LeapFrog_Term_free/HIPDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/HID/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_LeapFrog_Term_free/HIDDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/PHE/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_LeapFrog_Term_free/PHEDv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/TY3/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_LeapFrog_Term_free/TY3Dv_econ_300_100ps_1-10.econ.rmsd.av" , "/home/yamamori/calspa/TAMD_econ/dipep_2014-09-14_ABAMD_debugged/TRP/econ_wdi_NVE_AMBER_TermOn_tune_dihedafclust_and_min_wotuning_LeapFrog_Term_free/TRPDv_econ_300_100ps_1-10.econ.rmsd.av"  )
 econCMD <- c(	"/home/yamamori/calspa/refcalc-s-ivy3/refcalc/dipep/GLY/econ_wdi_NVE_100ps/GLYDv_econ_300_100ps_1-10_2ntc.econ.rmsd.av" , "/home/yamamori/calspa/refcalc-s-ivy3/refcalc/dipep/ALA/econ_wdi_NVE_100ps/ALADv_econ_300_100ps_1-10_2ntc.econ.rmsd.av" , "/home/yamamori/calspa/refcalc-s-ivy3/refcalc/dipep/VAL/econ_wdi_NVE_100ps/VALDv_econ_300_100ps_1-10_2ntc.econ.rmsd.av" , "/home/yamamori/calspa/refcalc-s-ivy3/refcalc/dipep/LEU/econ_wdi_NVE_100ps/LEUDv_econ_300_100ps_1-10_2ntc.econ.rmsd.av" , "/home/yamamori/calspa/refcalc-s-ivy3/refcalc/dipep/ILE/econ_wdi_NVE_100ps/ILEDv_econ_300_100ps_1-10_2ntc.econ.rmsd.av" , "/home/yamamori/calspa/refcalc-s-ivy3/refcalc/dipep/PRO/econ_wdi_NVE_100ps/PRODv_econ_300_100ps_1-10_2ntc.econ.rmsd.av" , "/home/yamamori/calspa/refcalc-s-ivy3/refcalc/dipep/MET/econ_wdi_NVE_100ps/METDv_econ_300_100ps_1-10_2ntc.econ.rmsd.av" , "/home/yamamori/calspa/refcalc-s-ivy3/refcalc/dipep/CYS/econ_wdi_NVE_100ps/CYSDv_econ_300_100ps_1-10_2ntc.econ.rmsd.av" , "/home/yamamori/calspa/refcalc-s-ivy3/refcalc/dipep/CYX/econ_wdi_NVE_100ps/CYXDv_econ_300_100ps_1-10_2ntc.econ.rmsd.av" , "/home/yamamori/calspa/refcalc-s-ivy3/refcalc/dipep/CYM/econ_wdi_NVE_100ps/CYMDv_econ_300_100ps_1-10_2ntc.econ.rmsd.av" , "/home/yamamori/calspa/refcalc-s-ivy3/refcalc/dipep/SER/econ_wdi_NVE_100ps/SERDv_econ_300_100ps_1-10_2ntc.econ.rmsd.av" , "/home/yamamori/calspa/refcalc-s-ivy3/refcalc/dipep/THR/econ_wdi_NVE_100ps/THRDv_econ_300_100ps_1-10_2ntc.econ.rmsd.av" , "/home/yamamori/calspa/refcalc-s-ivy3/refcalc/dipep/ASN/econ_wdi_NVE_100ps/ASNDv_econ_300_100ps_1-10_2ntc.econ.rmsd.av" , "/home/yamamori/calspa/refcalc-s-ivy3/refcalc/dipep/ASH/econ_wdi_NVE_100ps/ASHDv_econ_300_100ps_1-10_2ntc.econ.rmsd.av" , "/home/yamamori/calspa/refcalc-s-ivy3/refcalc/dipep/GLN/econ_wdi_NVE_100ps/GLNDv_econ_300_100ps_1-10_2ntc.econ.rmsd.av" , "/home/yamamori/calspa/refcalc-s-ivy3/refcalc/dipep/GLH/econ_wdi_NVE_100ps/GLHDv_econ_300_100ps_1-10_2ntc.econ.rmsd.av" , "/home/yamamori/calspa/refcalc-s-ivy3/refcalc/dipep/ASP/econ_wdi_NVE_100ps/ASPDv_econ_300_100ps_1-10_2ntc.econ.rmsd.av" , "/home/yamamori/calspa/refcalc-s-ivy3/refcalc/dipep/GLU/econ_wdi_NVE_100ps/GLUDv_econ_300_100ps_1-10_2ntc.econ.rmsd.av" , "/home/yamamori/calspa/refcalc-s-ivy3/refcalc/dipep/LYS/econ_wdi_NVE_100ps/LYSDv_econ_300_100ps_1-10_2ntc.econ.rmsd.av" , "/home/yamamori/calspa/refcalc-s-ivy3/refcalc/dipep/LYN/econ_wdi_NVE_100ps/LYNDv_econ_300_100ps_1-10_2ntc.econ.rmsd.av" , "/home/yamamori/calspa/refcalc-s-ivy3/refcalc/dipep/ARG/econ_wdi_NVE_100ps/ARGDv_econ_300_100ps_1-10_2ntc.econ.rmsd.av" , "/home/yamamori/calspa/refcalc-s-ivy3/refcalc/dipep/HIE/econ_wdi_NVE_100ps/HIEDv_econ_300_100ps_1-10_2ntc.econ.rmsd.av" , "/home/yamamori/calspa/refcalc-s-ivy3/refcalc/dipep/HIP/econ_wdi_NVE_100ps/HIPDv_econ_300_100ps_1-10_2ntc.econ.rmsd.av" , "/home/yamamori/calspa/refcalc-s-ivy3/refcalc/dipep/HID/econ_wdi_NVE_100ps/HIDDv_econ_300_100ps_1-10_2ntc.econ.rmsd.av" , "/home/yamamori/calspa/refcalc-s-ivy3/refcalc/dipep/PHE/econ_wdi_NVE_100ps/PHEDv_econ_300_100ps_1-10_2ntc.econ.rmsd.av" , "/home/yamamori/calspa/refcalc-s-ivy3/refcalc/dipep/TY3/econ_wdi_NVE_100ps/TY3Dv_econ_300_100ps_1-10_2ntc.econ.rmsd.av" , "/home/yamamori/calspa/refcalc-s-ivy3/refcalc/dipep/TRP/econ_wdi_NVE_100ps/TRPDv_econ_300_100ps_1-10_2ntc.econ.rmsd.av"  )

 file.name=paste("/home/yamamori/Report/2015-02/eps/fig_econ_vdt_100ps_wdi_NVE_ABAMD_AMBER_TermOn_sander_afclust_and_minimize_all_PC6vsLFvsCMD_2015-02-02.eps",sep='')
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
 iro    <-c( "red", "blue", "black" )
 file.names<-NULL

 for (i in 1:num.res) {
  par(mar=c(0.0,0.0,0.0,0.0))

  file.names[1]=econPC6[i]
  file.names[2]=econmLF[i]
  file.names[3]=econCMD[i]

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
  text(8,-7,RES[i],cex=1.75)
  box(lwd=2.0)

  if ( i >= 19 && i < 27 ) {
    xrange.axis <- c(0,15,15)
    axis(1,xaxp=xrange.axis,lwd=2.0,cex=0.75)
  }
  else if ( i == 27 ) {
    xrange.axis <- c(0,16,16)
    axis(1,xaxp=xrange.axis,lwd=2.0,cex=0.75)
  }

  if ( i == 1) {
    yrange.axis=c( -7,2,9 )
    axis(2,yaxp=yrange.axis,lwd=2.0,cex=0.75)
  }
  else if ( i == 10) {
    yrange.axis=c( -7,2,9 )
    axis(2,yaxp=yrange.axis,lwd=2.0,cex=0.75)
  }
  else if ( i == 19) {
    yrange.axis=c( -8,2,10 )
    axis(2,yaxp=yrange.axis,lwd=2.0,cex=0.75)
  }
    
}

mtext(outer=T,label.x,side=1,line=4.0,cex=1.0)
mtext(outer=T,label.y,side=2,line=3.0,cex=1.0)

