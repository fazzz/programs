fact.x <- 1
fact.y <- 1
fact.p <- 1

dirout <- "/home/yamamori/papers/CG-FG_TACCM_REMD"

title <- NULL

label.x <- " "
label.y <- " "

name.out <- paste(dirout,"/eps/","fig7_a-f_comp_MuSTARMD-REMD-REUS-TAMD-CMD_MetEnk_2Dpmf_2013-09-02",sep='')

pdatanamedummy <- paste(dirout,"/fig/points_dummy.txt",sep='')
pdataname <- paste(dirout,"/fig/points.txt",sep='')
ldataname <- paste(dirout,"/fig/lines.txt",sep='')
ldataname2 <- paste(dirout,"/fig/lines2.txt",sep='')
ldataname3 <- paste(dirout,"/fig/lines3.txt",sep='')
   
level <- seq(0,20,2)

xrange <- c(0,0.5,5)
xrange.axis <- c(0,0.5,5)
yrange <- c(0,0.5,5)
yrange.axis <- c(0.0,0.5,5)

file.name <- paste(name.out,'.eps',sep='')
postscript(file.name,width=5.3,height=6.0,horizontal=FALSE,onefile=FALSE,paper="special")

par(oma = c(0,0,0,0) )

source("~/papers/CG-FG_TACCM_REMD/FillConMap_wrange_worwoaxis.R")
source("~/papers/CG-FG_TACCM_REMD/FillConBar.R")

source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMap_wrange_worwoaxis_wpoints.R")
source("~/papers/CG-FG_TACCM_REMD/fel_wFillConBar_wmai.R")

source("~/papers/CG-FG_TACCM_REMD/fel_wFillConMapbox_wpoints_2.R")


nf <- layout(matrix(c(1,2,7,3,4,7,5,6,7),3,3,byrow=TRUE),c(22.5,16.5,8),c(20,15,23))

par(cex.axis=1.44)
par(cex.lab=1.44)

#name <- "~/calspa/TACCM_CGAAREMD/MetEnk/e_CG-FG_CMD_NH_2013-08-31/1FG2CG/tau=1.0/mZ=100.00/ep=0.01/cutoff=4.7/TZ=800/TCG=400/2CG_2013-09-02_KAA=1000_KCG=1000-500_2/freq=1/fb=1.0/fa=1.0/ft=0.2/fc=1.0/fn=1.0/pmf/pmf_TAA=300_TCG1_400_TCG2_400_TZ_800_0.03_1000_0_0_AA_bo10000ps_2013-08-24_2013-08-31"

name <- "~/calspa/TACCM_CGAAREMD/MetEnk/e_CG-FG_CMD_NH_2013-08-31/1FG2CG/tau=1.0/mZ=1000.00/ep=0.01/cutoff=4.7/TZ=1000/TCG=400/2CG_2013-09-05_KAA=1250_KCG=1000-500_2/freq=1/fb=1.0/fa=1.0/ft=0.2/fc=1.0/fn=1.0/pmf/pmf_TAA=300_TCG1_400_TCG2_400_TZ_1000_0.015_1250_0_0_AA_bo10000ps_2013-08-24_2013-08-31"

cat(name)

felwFillConMapwrangeworwoaxiswpoints(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,mai=c(0.0,0.8,0.5,0.0),xaflag="f",pdata=pdatanamedummy,iro=c("white","white","white"))

text(0.4,0.4,"(a)",cex=1.6)


name <- "~/calspa/refcalc/REMD/MetEnk/s_RE_v_2012-07-25_SBL1/f300t600_80ns/nEX10000/frq1/pmf/pmf_MetEnkv_0.025_0.5"

xrange <- c(0,0.5,5)
xrange.axis <- c(0.1,0.5,4)
yrange <- c(0,0.5,5)
yrange.axis <- c(0,0.4,4)

cat(name)

label.y <- ""

felwFillConMapwrangeworwoaxiswpoints(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,mai=c(0.0,0.0,0.5,0.1),xaflag="f",yaflag="f",pdata=pdatanamedummy,iro=c("white","white","white"))

text(0.4,0.4,"(b)",cex=1.6)


name <- "~/calspa/refcalc/REUS/MetEnk/s_REUSVAC_2013-08-30_ff99SB/UmbAt_Nbin_8_K_0.4_2/nEX_1000/freq_10ps/pmf/pmf_MetEnkv"

cat(name)

felwFillConMapwrangeworwoaxiswpoints(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,mai=c(0.0,0.8,0.0,0.0),xaflag="f",pdata=pdatanamedummy,iro=c("white","white","white"))

text(0.4,0.4,"(c)",cex=1.6)



name <- "~/calspa/TACCM/MetEnk/e_TACCM_NH_2013-08-28_99SB_80ns/tau=1.0/TB=600/KZ=1000/mZ=100.00/pmf/pmf_T=300.Zhist_0.01"

cat(name)

felwFillConMapwrangeworwoaxiswpoints(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,mai=c(0.0,0.0,0.0,0.1),xaflag="f",yaflag="f",pdata=pdatanamedummy,iro=c("white","white","white"))

text(0.4,0.4,"(d)",cex=1.6)


name <- "~/calspa/refcalc/CMD/MetEnk/s_CVAC_2013-03-13_LM1/anl/pmf_MetEnk_T=300_fLM1_0.01_0.0-0.5.txt"

cat(name)

felwFillConMapwrangeworwoaxiswpoints(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,mai=c(0.8,0.8,0.0,0.0),pdata=pdatanamedummy,iro=c("white","white","black"))

text(0.4,0.4,"(e)",cex=1.6)


name <- "~/calspa/refcalc/CMD/MetEnk/s_CVAC_2013-03-13_LM2/anl/pmf_MetEnk_T=300_fLM2_0.01_0.0-0.5.txt"

cat(name)

felwFillConMapwrangeworwoaxiswpoints(name,label.x=label.x,label.y=label.y,level=level,title=title,kt="F/kBT",xrange=xrange,yrange=yrange,xrange.axis=xrange.axis,yrange.axis=yrange.axis,mai=c(0.8,0.0,0.0,0.1),yaflag="f",pdata=pdatanamedummy,iro=c("white","white","black"))

text(0.4,0.4,"(f)",cex=1.6)


felwFillConBarwmai(name,label.x=label.x,label.y=label.y,level=level,mai=c(0.8,0.5,0.5,0.2))

