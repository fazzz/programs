source("~/Rspa/ini.R")
source("~/Rspa/set_spe.R")
source("~/Rspa/set_res.R")
source("~/Rspa/plmGeneral_wsrange.R")

hutosa <- 1

leg.pos="topright"
is.leg <- 0

xrange <- c(0,600,6)
yrange <- c(0,2.0,2)

name.res <- c( "ASN", "GLN", "LEU","SER", "LYS","ILE","ALA", "CYS", "MET", "THR", "VAL", "ARG", "ASP", "GLY", "HIS","PHE","PRO", "THR","TYR","TRP" )

num.res <- length(name.res)

tp <- c("30")
num.tp <- length(tp)

dirbase.nma <-"/misc/ivy3/yamamori/ajisai/TAMD/TAMD/spectrum/nmadas/"
dirbase.tamd <- "/misc/ivy3/yamamori/ajisai/TAMD/spectrum/tamd_ECEPP/testspace_debug_af_optq/"

file.name=paste("~/defense/eps/ABA/fig_spe_comp_20dipeptides_TAMDvsNMADAS_2013-07-09.eps",sep='')
postscript(file.name,width=8,height=6,horizontal=FALSE,onefile=FALSE,paper="special")

sen=c(1,2)
iro<-c(3,2)
id.ys<-c(2,2)
tenshu<-c(1,1)

par(mfrow=c(4,5))
par(mar=c(0.0,0.0,0.0,0.0))
par(oma=c(3.0,6.0,2.0,2.0))

for (i in 1:num.res) {
  file.names.nma=paste(dirbase.nma,name.res[i],'_3_13_af_sim_2_4/',name.res[i],'.freq_hist',sep='')
  file.names.tamd=paste(dirbase.tamd,name.res[i],'_spe_step_equ_2_4/',tp[1],'_0.1/spectrum/spe_',name.res[i],'Dv_c_a.txt',sep='')

  file.names=c(file.names.tamd,file.names.nma)

  plmGeneralwsrange(data.names=file.names,
                    id.ys=id.ys,
                    label.size=0.5,axis.size=2.0,
                    iro=iro,axis.ft="F",
                    xrange=xrange,yrange=yrange,
                    is.sen=sen,width=0.1,is.header="F")
  text(400,0.5,name.res[i],cex=1.5)
####################################################
##   if (i==1) {                                  ##
##     legend(200,1.9,legend=c("tamd"),bty="n",   ##
##            lty=c(1),cex=1.5,y.intersp=2,       ##
##            col=iro[1],ncol=leg.col)            ##
##     legend(200,1.7,legend=c("nmadas"),bty="n", ##
##            lty=c(1),cex=1.5,y.intersp=2,       ##
## #           cex=1.5,y.intersp=2,               ##
##            col=iro[2],ncol=leg.col)            ##
## #           ,pch=c(1))                         ##
##   }                                            ##
####################################################
  box(lwd=2.0)

  if (i==16 || i==17 || i==18 || i==19 || i==20) {
    xrange.axis <- c(0,500,5)
    axis(1,xaxp=xrange.axis,lwd=2.0)
  }

  if (i==20) {
    xrange.axis <- c(0,600,6)
    axis(1,xaxp=xrange.axis,lwd=2.0)
  }

  if (i==1 ) {
    yrange.axis <- c(0,2,2)
    axis(2,yaxp=yrange,lwd=2.0)
  }
  
  if (i==1 || i==6 || i==11 || i==16 ) {
    yrange.axis <- c(0,1,1)
    axis(2,yaxp=yrange,lwd=2.0)
  }

}

