source("~/Rspa/ini.R")
source("~/Rspa/set_enecon.sh")
source("~/Rspa/set_res.R")
source("~/Rspa/plmGeneral_wsrange.R")

leg.pos="top"
is.leg <- 1

name.res <- c( "ASN", "GLN", "LYS", "SER", "ALA", "ARG", "CYS", "ILE", "LEU", "MET", "THR", "VAL" )
num.res <- length(name.res)

#id.ys.gen <- c(8,9,7,14,9,11,7,10,9,12,10,12,10,9,5,9,10,10,9,6)

id.ys.gen <- c(7)
temp <- "300"
id.ys.wogen <- c(3,4,5,6,8,9,10,11)
n.id.ys.wogen <- length(id.ys.wogen)
temp <- "300"

hutosa <- 1

label.x="fraction of V_n"

xrange <- c(0,1.1,11)
xrange.axis <- c(0.1,1.0,9)

yrange <- c(0,700,7)
yrange.axis <- c(0,700,7)

name.out=paste("~/thesis/TAMD/eps/fig_Chapt5_move_spectrum_by_dihed_tune_12_dipeptides.eps",sep='')
postscript(name.out,width=6.0,height=5.0,horizontal=FALSE,onefile=FALSE,paper="special")

par(mfrow=c(4,5))
par(oma=c(7.5,4.0,2.0,2.0))
par(mar=c(0.0,0.0,0.0,0.0))

dirbase.name <-paste("~/calspa/TAMD-s-ivy3/TAMD/mod_parm/dipep",sep='')
dirbase2.name <- paste("spe_AMBER_TermOn_wdtune_dihed_V_n_term_2012-02-21_pc6",sep='')

label.names<-NULL
file.names<-NULL
id.ys <- NULL

iro[1] <- "red"
iro[2] <- "blue"

for ( i in 1:num.res  ) {
  j<-1
  cat(id.ys)
  sen <-c(3)
  id.ys <- 1
  file.names <- paste(dirbase.name,"/",name.res[i],"/spe_AMBER_TermOn_wdtune_dihed_V_n_term_2012-02-21_pc6/spe_dtrj_",name.res[i],"Dv_T=300_tau=1_1fs_100ps_max_period_sum.txt",sep='')
  label.names <- name.res[i]

  plmGeneralwsrange(data.names=file.names,
                    label.names=label.names,
                    id.ys=id.ys,
                    label.size=0.5,axis.size=2.0,
                    iro=iro,axis.ft="F",is.header="F",
                    sdiro=iro,
                    xrange=xrange,yrange=yrange,
                    sdyrange=yrange,
                    is.sen=sen,width=10.0,
                    ,warrow="F")
  text(300,0.4,name.res[i])
  box(lwd=2.0)

}

label.x <- expression(paste("frequency (",cm^{-1},")"))
mtext(outer=T,label.x,side=1,line=3.0,cex=1.0)        
label.y <- expression(paste("spectrum (cm)"))
mtext(outer=T,label.y,side=2,line=2.0,cex=1.0)

