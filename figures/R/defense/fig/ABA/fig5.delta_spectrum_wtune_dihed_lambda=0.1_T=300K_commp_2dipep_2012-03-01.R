source("~/Rspa/ini.R")

source("~/Rspa/set_enecon.sh")
source("~/Rspa/set_res.R")
source("~/Rspa/plmGeneral_wsrange.R")

leg.pos="topright"
is.leg <- 0

temp <- "300"
temp2 <- "300"

xrange <- c(0,9,9)

lambda <- "0.4"

name.res <- c( "ASN","GLN")

numini=10

num.res <- length(name.res)

dirbase.name <- "/home/yamamori/calspa/TAMD/dipeptides/"
dirbase21.name <-paste("~/calspa/TAMD/mod_parm/dipep",sep='')
dirbase22.name <- paste("econ_wdi_NVE_AMBER_TermOn_wtune_2012-02-29_pc6",sep='')

sen<-NULL
id.ys<-NULL
ids.ys<-NULL
tenshu<-NULL

for ( k in 1:2  ) {
  sen[k]<-c(1)
  id.ys[k]<-c(2)
  ids.ys[k]<-c(3)
  tenshu[k]<-c(20)
}

par(mfrow=c(3,2))
par(oma=c(7.5,4.0,2.0,2.0))

file.names<-NULL
k2<-1
par(mar=c(0.0,0.0,0.0,0.0))

dirbase.name2 <-paste("~/calspa/TAMD/mod_parm/dipep",sep='')
dirbase2.name2 <- paste("spe_AMBER_TermOn_wdtune_dihed_V_n_term_2012-02-21_pc6",sep='')
xrange2 <- c(0,600,6)
xrange.axis2 <- c(0,600,6)
  
label.names2<-NULL
file.names2<-NULL
  
sen2 <- NULL
iro2 <- NULL
sen2[1] <- 1
sen2[2] <- 1
senshu[1]<-1
senshu[2]<-2
  
file.names2 <- NULL
file.name.max.period2 <- NULL
id.ys2 <- NULL

i2 <- i/2
  
file.names2 <- NULL
file.name.max.period2 <- NULL
id.ys2 <- NULL
j2 <- i

for (i in 1:num.res) {
  file.names2 <- paste(dirbase.name2,"/",name.res[i],"/",dirbase2.name2,"/n=1.0/anl/spe_dtrj_",name.res[i],"Dv_T=",temp2,"_tau=1_1fs_100ps.txt",sep='')

  if ( i==1 )
    id.ys2 <- c(7)
  if ( i==2 )
    id.ys2 <- c(8)

  iro2[1]<-"red"
  
  yrange2 <- c(0,0.50,1)
  yrange.axis2 <- c(-1,1,1)
  
  cat(id.ys2,"\n")
  plmGeneralwsrange(data.names=file.names2,
                    id.ys=id.ys2,
                    label.size=0.5,axis.size=2.0,
                    iro=iro2,axis.ft="F",is.header="T",
                    sdiro=iro2,
                    xrange=xrange2,yrange=yrange2,
                    sdyrange=yrange2,
                    is.sen=sen,width=10.0,
                    ,warrow="F")
  box(lwd=2.0)
}  

for (i in 1:num.res) {
  file.names2 <- paste(dirbase.name2,"/",name.res[i],"/",dirbase2.name2,"/n=",lambda,"/anl/spe_dtrj_",name.res[i],"Dv_T=",temp2,"_tau=1_1fs_100ps.txt",sep='')

  if ( i==1 )
    id.ys2 <- c(7)
  if ( i==2 )
    id.ys2 <- c(8)

  iro2[1]<-"blue"

  yrange2 <- c(0,0.50,1)
  yrange.axis2 <- c(-1,1,1)
  
  cat(id.ys2,"\n")
  plmGeneralwsrange(data.names=file.names2,
                    id.ys=id.ys2,
                    label.size=0.5,axis.size=2.0,
                    iro=iro2,axis.ft="F",is.header="T",
                    sdiro=iro2,
                    xrange=xrange2,yrange=yrange2,
                    sdyrange=yrange2,
                    is.sen=sen,width=10.0,
                    ,warrow="F")
  box(lwd=2.0)
}  

for (i in 1:num.res) {
  file.names2 <- paste(dirbase.name2,"/",name.res[i],"/",dirbase2.name2,"/n=",lambda,"/anl/spe_diff_dtrj_",name.res[i],"Dv_T=",temp2,"_tau=1_1fs_100ps.txt",sep='')
      
  if ( i==1 )
    id.ys2 <- c(7)
  if ( i==2 )
    id.ys2 <- c(8)

  iro2[1]<-"black"
  
  yrange2 <- c(-0.20,0.20,1)
  yrange.axis2 <- c(-1,1,1)

  cat(id.ys2,"\n")
  plmGeneralwsrange(data.names=file.names2,
                    id.ys=id.ys2,
                    label.size=0.5,axis.size=2.0,
                    iro=iro2,axis.ft="F",is.header="T",
                    sdiro=iro2,
                    xrange=xrange2,yrange=yrange2,
                    sdyrange=yrange2,
                    is.sen=sen,width=10.0,
                    ,warrow="F")
  box(lwd=2.0)
  axis(1,xaxp=xrange.axis2,lwd=2.0)
  axis(2,yaxp=yrange.axis2,lwd=2.0)
}  

name.out=paste("~/gakkai/seibutu_kanto/fig5.delta_spectrum_wtune_dihed_lambda=0.1_T=300K_commp_2dipep_2012-03-01",sep='')
OutTiff(name.out,w_val=600,h_val=380)
#OutTiff(name.out,w_val=850,h_val=480)
