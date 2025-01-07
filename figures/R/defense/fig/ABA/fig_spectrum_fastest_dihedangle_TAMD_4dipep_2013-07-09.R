source("~/Rspa/ini.R")

source("~/Rspa/set_enecon.sh")
source("~/Rspa/set_res.R")
source("~/Rspa/plmGeneral_wsrange.R")

leg.pos="topright"
is.leg <- 0

temp <- "300"

#name.res <- c( "ASN","GLN","ASP","GLY")
name.res <- c("ASP","GLY","ASN","GLN")

numini=10

num.res <- length(name.res)

file.name=paste("~/defense/eps/ABA/fig_spectrum_fastest_dihedangle_TAMD_4dipep_2013-07-09.eps",sep='')
postscript(file.name,width=3.0,height=4.8,horizontal=FALSE,onefile=FALSE,paper="special")

par(mfrow=c(2,1))
par(oma=c(7.5,4.0,2.0,2.0))

dirbase.name2 <-paste("/misc/ivy3/yamamori/ajisai/TAMD/TAMD/mod_parm/dipep",sep='')
dirbase2.name2 <- paste("spe_AMBER_TermOn_wdtune_dihed_V_n_term_2012-02-21_pc6",sep='')

xrange <- c(0,600,6)
xrange.axis <- c(0,600,6)

yrange <- c(0,0.5,1)
yrange.axis <- c(-1,1,1)

label.names2<-NULL
file.names2<-NULL

hutosa <- 1

sen2 <- NULL
iro2 <- NULL
sen2[1] <- 1
sen2[2] <- 1
senshu[1]<-1
senshu[2]<-2
iro2[1]<-"black"
iro2[2]<-"black"

file.names2 <- NULL
file.name.max.period2 <- NULL
id.ys2 <- NULL

for (i in 1:2) {
  par(mar=c(0.0,1.0,0.0,0.0))
  i2 <- i/2

  file.names2 <- NULL
  file.name.max.period2 <- NULL
  id.ys2 <- NULL
  j2 <- i
  k2<-1
  for ( j in 1:2  ) {
    for ( j2 in 1:2  ) {
      file.names2[j2] <- paste(dirbase.name2,"/",name.res[k2],"/",dirbase2.name2,"/n=1.0/anl/spe_dtrj_",name.res[k2],"Dv_T=",temp,"_tau=1_1fs_100ps.txt",sep='')
      file.name.max.period2[j2] <- paste(dirbase.name2,"/",name.res[k2],"/",dirbase2.name2,"/n=1.0/anl/spe_dtrj_",name.res[k2],"Dv_T=",temp,"_tau=1_1fs_100ps_max_period.txt",sep='')

      b2 <- read.table(file.name.max.period2[j2])
      id.ys2[j2] <- b2[,1]
    k2 <- k2+1
    }
  
    plmGeneralwsrange(data.names=file.names2,
                      id.ys=id.ys2,
                      label.size=0.5,axis.size=2.0, 
                      iro=iro2,axis.ft="F",is.header="T",
                      sdiro=iro2,
                      xrange=xrange,yrange=yrange,
                      sdyrange=yrange2,
                      is.sen=sen,width=10.0,
                      warrow="F")           
    box(lwd=2.0)                             
    
    axis(2,yaxp=yrange.axis,lwd=2.0)
  }
  if (i==2)
    axis(1,xaxp=xrange.axis,lwd=2.0)
}

